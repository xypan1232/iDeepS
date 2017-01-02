#include "NSPDK_FeatureGenerator.h"

//----------------------------------------------------------------------------------------------------------------------------
void DebugClass::Clear() {
	mHashToFeatureMap.clear();
	mHashToGraphMap.clear();
}

void DebugClass::Output(ostream& out) const {
	OutputFeatureEncoding(out);
}

void DebugClass::OutputFeatureEncoding(ostream& out) const {
	out << "#Feature encodings: [" << mHashToFeatureMap.size() << "]" << endl;
	for (map<unsigned, string>::const_iterator it = mHashToFeatureMap.begin(); it != mHashToFeatureMap.end(); ++it)
		out << it->first << " -> " << it->second << endl;
}

void DebugClass::StoreFeatureCodeToFeatureInfo(unsigned aFeatureCode, vector<unsigned>& aDetailsList) {
	if (mHashToFeatureMap.count(aFeatureCode) == 0) {
		string feature_information = "r:" + stream_cast<string>(aDetailsList[0]) + " d:" + stream_cast<string>(aDetailsList[1]);
		feature_information += " [g: " + mHashToGraphMap[aDetailsList[2]] + "]";
		feature_information += " [g: " + mHashToGraphMap[aDetailsList[3]] + "]";
		mHashToFeatureMap.insert(make_pair(aFeatureCode, feature_information));
	}
}

void DebugClass::SerializedRootedGraphCanonicalFormEncoding(unsigned aDiscreteEncoding, int aRootVertexIndex, const GraphClass& aG, int aRadius) {
	string value = "";
	if (aRadius > 0) {
		//extract set of vertices in the ball of radius aMaxDepth
		set<unsigned> ball;
		for (int r = 0; r <= aRadius; r++) {
			vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(aRootVertexIndex, r);
			ball.insert(dest_id_list.begin(), dest_id_list.end());
		}
		if (ball.size() == 0)
			throw std::logic_error("ERROR99: Something went wrong in SerializedRootedGraphCanonicalFormEncoding: cannot generate features over an empty neighbourhood graph!");
		//induce the subgraph from the ball and return the new index for the root vertex
		GraphClass gal;
		unsigned root = aG.GetVertexInducedRootedSubGraph(ball, aRootVertexIndex, gal);
		string root_label = gal.GetVertexLabel(root);
		root_label += "*";
		gal.SetVertexLabel(root, root_label);
		value = gal.Serialize();
	} else
		value = "*" + aG.GetVertexLabelConcatenated(aRootVertexIndex);
	mHashToGraphMap.insert(make_pair(aDiscreteEncoding, value));
}

//----------------------------------------------------------------------------------------------------------------------------------------------
NSPDK_FeatureGenerator::NSPDK_FeatureGenerator(const std::string& id) :
		FeatureGenerator(id), FlagsServiceClient(id) {
	mRadius = 0;
	mDistance = 0;
	mMatchType = "hard";
	mHashBitSize = (unsigned) (numeric_limits<unsigned>::digits - 1);
	mHashBitMask = numeric_limits<unsigned>::max() >> 1;
	mMinKernel = false;
	mNormalization = true;
	mUseRealVectorInformation = false;
	mNumHashFunctionsMinHash = 10;
	mDebugVerbosity = 0;
	mVertexDegreeThreshold = 10;
	new_flag(&mVertexDegreeThreshold, "vertex_degree_threshold", "(unsigned)\nThreshold vertex degree above which features for a specific vertex are not generated");
	new_flag(&mRadius, "radius", "(unsigned)\nMax radius of kernel neighborhoods");
	new_flag(&mDistance, "distance", "(unsigned)\nMax distance between pairs of neighborhoods");
	new_flag(&mMatchType, "match_type", "(string)\nHow to match neighborhoods: soft, hard");
	new_flag(&mNormalization, "normalization", "(bool)\nNormalize feature vectors");
	new_flag(&mMinKernel, "min_kernel", "(bool)\nApply the min-kernel on top of the generated features");
	new_flag(&mNumHashFunctionsMinHash, "num_hash_functions_minhash", "(unsigned)\nNumber of hash funcions in the minhash signature");

	new_flag(&mUseRealVectorInformation, "use_real_vector_information", "(bool)\nUse information encoded as a real vector per vertex");
	new_flag(&mHashBitSize, "hash_bit_size", "(unsigned)\nNumber of bits for hash values"); //FIXME: since parameter variables are accessed directly the bit size is useless as setting it cannot trigger automatically the computation of the bit mask; the only solution is to set directly the bit mask itself
	mHashBitMask = (2 << (mHashBitSize - 1)) - 1;
	new_flag(&mHashBitMask, "hash_bit_mask", "(unsigned)\nMask for hash values");
	new_flag(&mDebugVerbosity, "verbosity", "(unsigned)\nNumber for debug verbosity level");
}

void NSPDK_FeatureGenerator::OutputParameters(ostream& out) const {
	out << "Radius: " << mRadius << endl;
	out << "Distance: " << mDistance << endl;
	out << "Match_Type: " << mMatchType << endl;
	out << "Hash_Bit_Size: " << mHashBitSize << endl;
	out << "Hash_Bit_mask: " << mHashBitMask << endl;
	out << "Min_Kernel: " << mMinKernel << endl;
	out << "Normalization: " << mNormalization << endl;
	out << "Use_Real_Vector_Information: " << mUseRealVectorInformation << endl;
	out << "Num_hash_functions_minhash: " << mNumHashFunctionsMinHash << endl;
	out << "Vertex Degree Threshold: " << mVertexDegreeThreshold << endl;
	out << "Debug_Verbosity: " << mDebugVerbosity << endl;
}

void NSPDK_FeatureGenerator::OutputFeatureMap(ostream& out) const {
	mDebugInfo.OutputFeatureEncoding(out);
}

unsigned NSPDK_FeatureGenerator::GenerateGraphHashCode(const GraphClass& aG, const vector<unsigned>& aFirstEndpointList) {
	SVector x;
	generate_feature_vector(aG, x, aFirstEndpointList);
	unsigned code = GenerateVectorHashCode(x);
	return code;
}

void NSPDK_FeatureGenerator::InitFeatureCache(const GraphClass& aG, unsigned aRadius) {
	mFeatureCache.clear();
	unsigned vertex_size = aG.VertexSize();
	for (unsigned i = 0; i <= aRadius; ++i) {
		mFeatureCache.push_back(vector<unsigned>(vertex_size, 0));
	}
}

void NSPDK_FeatureGenerator::generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList) {
	vector<unsigned> first_endpoint_list = aFirstEndpointList;
	GetFirstEndpoints(aG, first_endpoint_list);
	if (first_endpoint_list.size() == 0)
		throw std::logic_error("ERROR9: Something went wrong: cannot generate features over an empty set of first endpoints!");
	int horizon = max(mDistance, mRadius);
	aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
	if (aG.Check() == false)
		throw logic_error("ERROR10: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

	InitFeatureCache(aG, mRadius);
	for (unsigned r = 0; r <= mRadius; r++) {
		for (unsigned d = 0; d <= mDistance; d++) {
			SVector z;
			for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
				unsigned src_id = first_endpoint_list[i];
				if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
					SVector zv;
					GenerateVertexFeatures(src_id, aG, r, d, zv);
					z.add(zv);
				}
			}
			if (mNormalization)
				z.normalize();
			x.add(z);
		}
	}

	if (mNormalization)
		x.normalize();
	if (mMinKernel)
		ConvertSparseVectorToMinFeatureVector(x);
	if (mDebugVerbosity > 0) {
		cout << x << endl;
		OutputFeatureMap(cout);
		aG.Output(cout);
	}
}

void NSPDK_FeatureGenerator::generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList) {
	vector<unsigned> first_endpoint_list = aFirstEndpointList;
	GetFirstEndpoints(aG, first_endpoint_list);
	if (first_endpoint_list.size() == 0)
		throw std::logic_error("ERROR6: Something went wrong: cannot generate features over an empty set of first endpoints!");
	int horizon = max(mDistance, mRadius);
	aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
	if (aG.Check() == false)
		throw logic_error("ERROR11: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

	InitFeatureCache(aG, mRadius);
	for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
		SVector z;
		unsigned src_id = first_endpoint_list[i];
		if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
			for (unsigned r = 0; r <= mRadius; r++) {
				for (unsigned d = 0; d <= mDistance; d++) {
					SVector zv;
					GenerateVertexFeatures(src_id, aG, r, d, zv);
					z.add(zv);
				}
			}
		}
		if (mMinKernel)
			ConvertSparseVectorToMinFeatureVector(z);
		x_list.push_back(z);
	}
	if (mDebugVerbosity > 0) {
		for (unsigned i = 0; i < x_list.size(); ++i)
			cout << i << " " << x_list[i] << endl;
		OutputFeatureMap(cout);
		aG.Output(cout);
	}
}

unsigned NSPDK_FeatureGenerator::GenerateVectorHashCode(SVector& x) {
	vector<pair<int, double> > vec = x.unpack();
	vector<unsigned> hash_vec;
	for (unsigned i = 0; i < vec.size(); ++i) {
		int key = vec[i].first;
		double val = vec[i].second;
		hash_vec.push_back((unsigned) key);
		hash_vec.push_back((unsigned) (val * 10000)); //Note: in general the value associated to a feature is a real number and it is therefore converted into an unsigned via scaling
	}
	unsigned code = HashFunc(hash_vec, mHashBitMask);
	return code;
}

unsigned NSPDK_FeatureGenerator::GenerateVertexHashCode(unsigned aSrcID, const GraphClass& aG) {
	SVector x;
	GenerateVertexFeatures(aSrcID, aG, x);
	unsigned code = GenerateVectorHashCode(x);
	return code;
}

void NSPDK_FeatureGenerator::GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, SVector& x) {
	for (unsigned r = 0; r <= mRadius; r++) {
		for (unsigned d = 0; d <= mDistance; d++) {
			SVector zv;
			GenerateVertexFeatures(aSrcID, aG, r, d, zv);
			x.add(zv);
		}
	}
}

void NSPDK_FeatureGenerator::GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x) {
	vector<unsigned> endpoint_list(4);
	endpoint_list[0] = aRadius;
	endpoint_list[1] = aDistance;

	unsigned src_code = GenerateVertexNeighbourhoodHashCode(aSrcID, aG, aRadius);
	vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(aSrcID, aDistance);
	for (unsigned dest_j = 0; dest_j < dest_id_list.size(); dest_j++) {
		unsigned dest_id = dest_id_list[dest_j];
		unsigned dest_code = 0;
		if (aG.GetVertexKernelPoint(dest_id) && aG.GetVertexAlive(dest_id)) { //proceed to extract features only if the *dest* vertex is a kernel point and is alive
			dest_code = GenerateVertexNeighbourhoodHashCode(dest_id, aG, aRadius);
			//impose canonical order for pair: i.e. A-B and B-A must generate the same feature
			if (src_code < dest_code) {
				endpoint_list[2] = src_code;
				endpoint_list[3] = dest_code;
			} else {
				endpoint_list[2] = dest_code;
				endpoint_list[3] = src_code;
			}
			unsigned code = HashFunc(endpoint_list, mHashBitMask);
			if (mDebugVerbosity > 0)
				mDebugInfo.StoreFeatureCodeToFeatureInfo(code, endpoint_list);
			SVector z;
			z.set(code, 1);
			x.add(z);
		}
	}
	if (mUseRealVectorInformation) {
		//get real vector per vertex
		vector<double> real_vector = aG.GetVertexNumericAttributeList(aSrcID);
		const unsigned num_hash_functions = mNumHashFunctionsMinHash;
		//get minhash of sparse feature vector per vertex
		vector<unsigned> signature = ComputeMinHashSignature(x, num_hash_functions);

		//convolute real vector with sparse feature vector per vertex
		//i.e. offset a copy of the real vector for each element of the minhash signature
		SVector z;
		vector<unsigned> convolution_feature_id(3);
		for (unsigned k = 0; k < signature.size(); k++) {
			convolution_feature_id[0] = k;
			vector<pair<int, double> > vec = x.unpack();
			for (unsigned i = 0; i < vec.size(); i++) {
				unsigned original_feature_id = vec[i].first;
				convolution_feature_id[1] = original_feature_id;
				for (unsigned j = 0; j < real_vector.size(); ++j) {
					unsigned original_real_vector_feature_id = j;
					double original_real_vector_value = real_vector[j];
					convolution_feature_id[2] = original_real_vector_feature_id;
					unsigned code = HashFunc(convolution_feature_id, mHashBitMask);
					SVector y;
					y.set(code, original_real_vector_value);
					z.add(y);
				}
			}
		}
		x = z;
	}
}

unsigned NSPDK_FeatureGenerator::GenerateVertexNeighbourhoodHashCode(unsigned aSrcID, const GraphClass& aG, unsigned aRadius) {
	unsigned src_code = 1;
	if (mFeatureCache[aRadius][aSrcID] == 0) {
		src_code = RootedGraphCanonicalFormEncoding(aSrcID, aG, aRadius);
		mFeatureCache[aRadius][aSrcID] = src_code;
	} else {
		src_code = mFeatureCache[aRadius][aSrcID];
	}
	return src_code;
}

void NSPDK_FeatureGenerator::ConvertSparseVectorToMinFeatureVector(SVector& x) {
	vector<pair<int, double> > vec = x.unpack();
	vector<unsigned> hash_vec(2, 0);
	SVector z;
	for (unsigned i = 0; i < vec.size(); ++i) {
		int key = vec[i].first;
		double val = vec[i].second;
		hash_vec[0] = (unsigned) (key);
		for (unsigned j = 0; j < val; ++j) {
			hash_vec[1] = j;
			unsigned code = HashFunc(hash_vec, mHashBitMask);
			z.set(code, 1); //NOTE: unresolved issues in case of collisions
		}
	}
	x = z;
}

unsigned NSPDK_FeatureGenerator::RootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG, int aRadius) {
	//NOTE:for efficiency reasons case radius=0 and radius=1 are treated as special cases
	unsigned discrete_encoding = 1;
	unsigned root_degree = aG.VertexAdjacentListSize(aRootVertexIndex);
	if (root_degree > mVertexDegreeThreshold)
		return discrete_encoding;

	if (aRadius == 0) { //return root label
		discrete_encoding = Radius0RootedGraphCanonicalFormEncoding(aRootVertexIndex, aG);
	} else if (aRadius == 1) { //return the sorted sequence of root's children
		discrete_encoding = Radius1RootedGraphCanonicalFormEncoding(aRootVertexIndex, aG);
	} else { //general case
		discrete_encoding = RadiusKRootedGraphCanonicalFormEncoding(aRootVertexIndex, aG, aRadius);
	}
	if (mDebugVerbosity > 0)
		mDebugInfo.SerializedRootedGraphCanonicalFormEncoding(discrete_encoding, aRootVertexIndex, aG, aRadius);
	return discrete_encoding;
}

unsigned NSPDK_FeatureGenerator::Radius0RootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG) {
	string encoding = aG.GetVertexLabelConcatenated(aRootVertexIndex);
	unsigned hash_subgraph_code = HashFunc(encoding);
	return hash_subgraph_code;
}

unsigned NSPDK_FeatureGenerator::Radius1RootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG) {
	unsigned hash_subgraph_code = 1;
	string encoding;
	encoding = aG.GetVertexLabelConcatenated(aRootVertexIndex) + ":";
	vector<pair<string, unsigned> > vertex_label_id_list;
	vector<unsigned> vertex_adjacency_list = aG.GetVertexAdjacentList(aRootVertexIndex);
	vector<unsigned> edge_adjacency_list = aG.GetEdgeAdjacentList(aRootVertexIndex);
	for (unsigned i = 0; i < vertex_adjacency_list.size(); ++i) {
		unsigned child_vertex_id = vertex_adjacency_list[i];
		unsigned child_edge_id = edge_adjacency_list[i];
		string child_label = aG.GetVertexLabelConcatenated(child_vertex_id) + "-" + aG.GetEdgeLabelConcatenated(child_edge_id);
		vertex_label_id_list.push_back(make_pair(child_label, child_vertex_id));
	}
	sort(vertex_label_id_list.begin(), vertex_label_id_list.end());
	if (vertex_label_id_list.size() > 0)
		encoding += vertex_label_id_list[0].first;
	for (unsigned i = 1; i < vertex_label_id_list.size(); i++)
		encoding += "." + vertex_label_id_list[i].first;
	hash_subgraph_code = HashFunc(encoding);
	return hash_subgraph_code;
}

unsigned NSPDK_FeatureGenerator::RadiusKRootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG, int aRadius) {
	unsigned hash_subgraph_code = 1;
	//extract set of vertices in the ball of radius aRadius
	set<unsigned> ball;
	for (int r = 0; r <= aRadius; r++) {
		vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(aRootVertexIndex, r);
		ball.insert(dest_id_list.begin(), dest_id_list.end());
	}

	//induce the subgraph from the ball and return the new index for the root vertex
	GraphClass gal;
	unsigned root = aG.GetVertexInducedRootedSubGraph(ball, aRootVertexIndex, gal);
	gal.ComputePairwiseDistanceInformation(gal.VertexSize());

	hash_subgraph_code = RootedGraphCanonicalFormEncoding(gal, root);
	return hash_subgraph_code;
}

void NSPDK_FeatureGenerator::CanonicalGraphVertexList(const GraphClass& aG, vector<unsigned>& oVertexList) {
	oVertexList.clear();
	GraphClass g(aG);
	g.ComputePairwiseDistanceInformation(g.VertexSize());
	for (unsigned i = 0; i < g.VertexSize(); ++i) {
		oVertexList.push_back(RootedGraphCanonicalFormEncoding(g, i));
	}
}

unsigned NSPDK_FeatureGenerator::MultiRootedGraphCanonicalFormEncoding(vector<int>& aRootVertexIndexList, const GraphClass& aG) {
	GraphClass g(aG); //make a copy of the graph
	g.ComputePairwiseDistanceInformation(g.VertexSize()); //compute complete breadth first visit

	//for all vertices extract the vertex's signature: sorted distance from all the other vertices + their vertex label + distance from labeled root vertices
	vector<unsigned> vertex_encoding_list;
	for (unsigned i = 0; i < g.VertexSize(); ++i) {
		vector<unsigned> vertex_encoding;
		vector<unsigned> vertex_distance_list;

		for (unsigned t = 0; t < aRootVertexIndexList.size(); ++t) {
			int dist = aG.PairwiseDistance(aRootVertexIndexList[t], i);
			string distance_label = stream_cast<string>(dist) + aG.GetVertexLabelConcatenated(aRootVertexIndexList[t]);
			unsigned hash_distance_label = HashFunc(distance_label);
			vertex_distance_list.push_back(hash_distance_label);
		}

		for (unsigned j = 0; j < g.VertexSize(); ++j) {
			int dist = g.PairwiseDistance(i, j);
			string distance_label = stream_cast<string>(dist) + g.GetVertexLabelConcatenated(j);
			unsigned hash_distance_label = HashFunc(distance_label);
			vertex_distance_list.push_back(hash_distance_label);
		}
		sort(vertex_distance_list.begin(), vertex_distance_list.end());
		unsigned hash_encoding = HashFunc(vertex_distance_list);
		vertex_encoding_list.push_back(hash_encoding);
	}
	//extract list of all edge's signatures in the induced graph: v-signature,u-signature,label(uv)
	vector<unsigned> edge_list;
	vector<unsigned> edge_encoding(3);
	for (unsigned u = 0; u < g.VertexSize(); ++u) {
		//get all edges of vertex u
		vector<unsigned> vertex_adjacency_list = g.GetVertexAdjacentList(u);
		vector<unsigned> edge_adjacency_list = g.GetEdgeAdjacentList(u);
		if (vertex_adjacency_list.size() == 0) { //NOTE: if a vertex is isolated use its encoding as edge encoding
			edge_list.push_back(vertex_encoding_list[u]);
		}
		for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {
			unsigned v = vertex_adjacency_list[j];

			if (vertex_encoding_list[u] < vertex_encoding_list[v]) {
				edge_encoding[0] = vertex_encoding_list[u];
				edge_encoding[1] = vertex_encoding_list[v];
			} else {
				edge_encoding[0] = vertex_encoding_list[v];
				edge_encoding[1] = vertex_encoding_list[u];
			}
			unsigned e = edge_adjacency_list[j];
			string edge_label = g.GetEdgeLabelConcatenated(e);
			unsigned hash_edge_label = HashFunc(edge_label);
			edge_encoding[2] = hash_edge_label;
			unsigned hash_edge_encoding = HashFunc(edge_encoding);
			edge_list.push_back(hash_edge_encoding);
		}
	}
	//the graph encoding is the sorted list of edge encodings
	sort(edge_list.begin(), edge_list.end());
	unsigned hash_subgraph_code = HashFunc(edge_list);
	return hash_subgraph_code;
}

unsigned NSPDK_FeatureGenerator::GraphCanonicalFormEncoding(const GraphClass& aG) {
	GraphClass g(aG); //make a copy of the graph
	g.ComputePairwiseDistanceInformation(g.VertexSize()); //compute complete breadth first visit

	//for all vertices extract the vertex's signature: sorted distance from all the other vertices + their vertex label
	vector<unsigned> vertex_encoding_list;
	for (unsigned i = 0; i < g.VertexSize(); ++i) {
		vector<unsigned> vertex_encoding;
		vector<unsigned> vertex_distance_list;
		for (unsigned j = 0; j < g.VertexSize(); ++j) {
			int dist = g.PairwiseDistance(i, j);
			string distance_label = stream_cast<string>(dist) + g.GetVertexLabelConcatenated(j);
			unsigned hash_distance_label = HashFunc(distance_label);
			vertex_distance_list.push_back(hash_distance_label);
		}
		sort(vertex_distance_list.begin(), vertex_distance_list.end());
		unsigned hash_encoding = HashFunc(vertex_distance_list);
		vertex_encoding_list.push_back(hash_encoding);
	}
	//extract list of all edge's signatures in the induced graph: v-signature,u-signature,label(uv)
	vector<unsigned> edge_list;
	vector<unsigned> edge_encoding(3);
	for (unsigned u = 0; u < g.VertexSize(); ++u) {
		//get all edges of vertex u
		vector<unsigned> vertex_adjacency_list = g.GetVertexAdjacentList(u);
		vector<unsigned> edge_adjacency_list = g.GetEdgeAdjacentList(u);
		if (vertex_adjacency_list.size() == 0) { //NOTE: if a vertex is isolated use its encoding as edge encoding
			edge_list.push_back(vertex_encoding_list[u]);
		}
		for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {
			unsigned v = vertex_adjacency_list[j];

			if (vertex_encoding_list[u] < vertex_encoding_list[v]) {
				edge_encoding[0] = vertex_encoding_list[u];
				edge_encoding[1] = vertex_encoding_list[v];
			} else {
				edge_encoding[0] = vertex_encoding_list[v];
				edge_encoding[1] = vertex_encoding_list[u];
			}
			unsigned e = edge_adjacency_list[j];
			string edge_label = g.GetEdgeLabelConcatenated(e);
			unsigned hash_edge_label = HashFunc(edge_label);
			edge_encoding[2] = hash_edge_label;
			unsigned hash_edge_encoding = HashFunc(edge_encoding);
			edge_list.push_back(hash_edge_encoding);
		}
	}
	//the graph encoding is the sorted list of edge encodings
	sort(edge_list.begin(), edge_list.end());
	unsigned hash_subgraph_code = HashFunc(edge_list);
	return hash_subgraph_code;
}

unsigned NSPDK_FeatureGenerator::RootedGraphCanonicalFormEncoding(const GraphClass& aG, unsigned aRootID) {
	//for all vertices extract the vertex's signature: distance from root-sorted distance from all the other vertices + their vertex label
	vector<unsigned> vertex_encoding_list;
	for (unsigned i = 0; i < aG.VertexSize(); ++i) {
		vector<unsigned> vertex_distance_list;
		int dist = aG.PairwiseDistance(aRootID, i);
		string distance_label = stream_cast<string>(dist) + aG.GetVertexLabelConcatenated(aRootID);
		unsigned hash_distance_label = HashFunc(distance_label);
		vertex_distance_list.push_back(hash_distance_label);

		for (unsigned j = 0; j < aG.VertexSize(); ++j) {
			int dist = aG.PairwiseDistance(i, j);
			string distance_label = stream_cast<string>(dist) + aG.GetVertexLabelConcatenated(j);
			unsigned hash_distance_label = HashFunc(distance_label);
			vertex_distance_list.push_back(hash_distance_label);
		}
		sort(vertex_distance_list.begin(), vertex_distance_list.end());
		unsigned hash_encoding = HashFunc(vertex_distance_list);
		vertex_encoding_list.push_back(hash_encoding);
	}
	//extract list of all edge's signatures in the induced graph: v-signature,u-signature,label(uv)
	vector<unsigned> edge_list;
	vector<unsigned> edge_encoding(3);
	for (unsigned u = 0; u < aG.VertexSize(); ++u) {
		//get all edges of vertex u
		vector<unsigned> vertex_adjacency_list = aG.GetVertexAdjacentList(u);
		vector<unsigned> edge_adjacency_list = aG.GetEdgeAdjacentList(u);
		if (vertex_adjacency_list.size() == 0) { //NOTE: if a vertex is isolated use its encoding as edge encoding
			edge_list.push_back(vertex_encoding_list[u]);
		}
		for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {
			unsigned v = vertex_adjacency_list[j];

			if (vertex_encoding_list[u] < vertex_encoding_list[v]) {
				edge_encoding[0] = vertex_encoding_list[u];
				edge_encoding[1] = vertex_encoding_list[v];
			} else {
				edge_encoding[0] = vertex_encoding_list[v];
				edge_encoding[1] = vertex_encoding_list[u];
			}
			unsigned e = edge_adjacency_list[j];
			string edge_label = aG.GetEdgeLabelConcatenated(e);
			unsigned hash_edge_label = HashFunc(edge_label);
			edge_encoding[2] = hash_edge_label;
			unsigned hash_edge_encoding = HashFunc(edge_encoding);
			edge_list.push_back(hash_edge_encoding);
		}
	}
	//the graph encoding is the sorted list of edge encodings
	sort(edge_list.begin(), edge_list.end());
	unsigned hash_subgraph_code = HashFunc(edge_list);
	return hash_subgraph_code;
}

void NSPDK_FeatureGenerator::GetFirstEndpoints(const GraphClass& aG, vector<unsigned>& oFirstEndpointList) const {
	//insert additional vertices
	if (oFirstEndpointList.size() == 0) //if oFirstEndpointList is empty then fill it with all vertices that are viewpoints, otherwise do nothing i.e. use the given list
		for (unsigned i = 0; i < aG.VertexSize(); ++i)
			if (aG.GetVertexViewPoint(i) || aG.GetVertexAbstraction(i))
				oFirstEndpointList.push_back(i);
	if (mDebugVerbosity > 0) {
		cout << "First endpoint id list [" << oFirstEndpointList.size() << "]:" << endl;
		for (unsigned i = 0; i < oFirstEndpointList.size(); i++)
			cout << oFirstEndpointList[i] << " ";
		cout << endl;
	}
}

inline
unsigned NSPDK_FeatureGenerator::HashFunc(const string& aString, unsigned aBitMask) { //NOTE: extract the least significant bits from the hash
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aString.length(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aString[i] * (hash >> 3)) : (~(((hash << 11) + aString[i]) ^ (hash >> 5)));
	}
	return hash & aBitMask;
}

inline
unsigned NSPDK_FeatureGenerator::HashFunc(const vector<unsigned>& aList, unsigned aBitMask) {
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aList.size(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aList[i] * (hash >> 3)) : (~(((hash << 11) + aList[i]) ^ (hash >> 5)));
	}
	return hash & aBitMask;
}

//----------------------------------------------------------------------------------------------------------------------------------------------
SK_FeatureGenerator::SK_FeatureGenerator(const std::string& id) :
		NSPDK_FeatureGenerator(id) {
}
void SK_FeatureGenerator::generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList){
	InitShellCache(aG,mRadius);
	NSPDK_FeatureGenerator::generate_feature_vector(aG,x,aFirstEndpointList);
}

void SK_FeatureGenerator::generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList){
	InitShellCache(aG,mRadius);
	NSPDK_FeatureGenerator::generate_vertex_feature_vector(aG,x_list,aFirstEndpointList);
}

void SK_FeatureGenerator::GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x) {
	vector<unsigned> endpoint_list(5);
	endpoint_list[0] = aRadius;
	endpoint_list[1] = aDistance;

	for (unsigned inner_radius = 0; inner_radius <= aRadius; ++inner_radius) {
		endpoint_list[4] = inner_radius;
		unsigned src_code = GenerateOuterRadiusKInnerRadiusMVertexNeighbourhoodHashCode(aSrcID, aG, aRadius, inner_radius);
		vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(aSrcID, aDistance);
		for (unsigned dest_j = 0; dest_j < dest_id_list.size(); dest_j++) {
			unsigned dest_id = dest_id_list[dest_j];
			unsigned dest_code = 0;
			if (aG.GetVertexKernelPoint(dest_id) && aG.GetVertexAlive(dest_id)) { //proceed to extract features only if the *dest* vertex is a kernel point and is alive
				dest_code = GenerateOuterRadiusKInnerRadiusMVertexNeighbourhoodHashCode(dest_id, aG, aRadius, inner_radius);
				//impose canonical order for pair: i.e. A-B and B-A must generate the same feature
				if (src_code < dest_code) {
					endpoint_list[2] = src_code;
					endpoint_list[3] = dest_code;
				} else {
					endpoint_list[2] = dest_code;
					endpoint_list[3] = src_code;
				}
				unsigned code = HashFunc(endpoint_list, mHashBitMask);
				if (mDebugVerbosity > 0)
					mDebugInfo.StoreFeatureCodeToFeatureInfo(code, endpoint_list);
				SVector z;
				z.set(code, 1);
				x.add(z);
			}
		}
	}
}

void SK_FeatureGenerator::InitShellCache(const GraphClass& aG, unsigned aRadius) {
	mShellCache.clear();
	unsigned vertex_size = aG.VertexSize();
	for (unsigned i = 0; i <= aRadius; ++i) {
		vector<vector<unsigned> > inner_vec;
		for (unsigned j = 0; j <= i; ++j)
			inner_vec.push_back(vector<unsigned>(vertex_size, 0));
		mShellCache.push_back(inner_vec);
	}
}

unsigned SK_FeatureGenerator::GenerateOuterRadiusKInnerRadiusMVertexNeighbourhoodHashCode(int aRootVertexIndex, const GraphClass& aG, int aOuterRadius, int aInnerRadius) {
	unsigned src_code = 1;
	if (mShellCache[aOuterRadius][aInnerRadius][aRootVertexIndex] == 0) {
		src_code = OuterRadiusKInnerRadiusMGraphCanonicalFormEncoding(aRootVertexIndex, aG, aOuterRadius, aInnerRadius);
		mShellCache[aOuterRadius][aInnerRadius][aRootVertexIndex] = src_code;
	} else {
		src_code = mShellCache[aOuterRadius][aInnerRadius][aRootVertexIndex];
	}
	return src_code;
}


unsigned SK_FeatureGenerator::OuterRadiusKInnerRadiusMGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG, int aOuterRadius, int aInnerRadius) {
	unsigned hash_subgraph_code = 1;
	//extract set of vertices in the ball of radius aRadius
	set<unsigned> ball;
	for (int r = aInnerRadius; r <= aOuterRadius; r++) {
		vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(aRootVertexIndex, r);
		ball.insert(dest_id_list.begin(), dest_id_list.end());
	}

	//induce the subgraph from the ball and return the new index for the root vertex
	GraphClass gal;
	vector<unsigned> original_vertex_map = aG.GetVertexInducedSubGraph(ball, gal);
	gal.ComputePairwiseDistanceInformation(gal.VertexSize());
	hash_subgraph_code = GraphCanonicalFormEncoding(gal);
	return hash_subgraph_code;
}

//----------------------------------------------------------------------------------------------------------------------------------------------
USPK_FeatureGenerator::USPK_FeatureGenerator(const std::string& id) :
		NSPDK_FeatureGenerator(id) {
}

void USPK_FeatureGenerator::GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x) {
	vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(aSrcID, aDistance);
	for (unsigned dest_j = 0; dest_j < dest_id_list.size(); dest_j++) {
		unsigned dest_id = dest_id_list[dest_j];
		{
			//extract graph at aSrc end
			set<unsigned> vertex_list = aG.GetUnionShortestPathsVertexIDList(aSrcID, dest_id, aDistance, aRadius);
			GraphClass gal;
			unsigned root = aG.GetVertexInducedRootedSubGraph(vertex_list, aSrcID, gal);
			gal.ComputePairwiseDistanceInformation(aRadius);
			unsigned code = RootedGraphCanonicalFormEncoding(gal, root);
			SVector z;
			z.set(code, 1);
			x.add(z);
		}
		{
			//extract graph at dest_id end
			set<unsigned> vertex_list = aG.GetUnionShortestPathsVertexIDList(dest_id, aSrcID, aDistance, aRadius);
			GraphClass gal;
			unsigned root = aG.GetVertexInducedRootedSubGraph(vertex_list, dest_id, gal);
			gal.ComputePairwiseDistanceInformation(aRadius);
			unsigned code = RootedGraphCanonicalFormEncoding(gal, root);
			SVector z;
			z.set(code, 1);
			x.add(z);
		}
	}
}

//----------------------------------------------------------------------------------------------------------------------------------------------
WDK_FeatureGenerator::WDK_FeatureGenerator(const std::string& id) :
		NSPDK_FeatureGenerator(id) {
	mLowerVertexDegreeThreshold = 1;
	new_flag(&mLowerVertexDegreeThreshold, "lower_vertex_degree_threshold", "(unsigned)\nThreshold vertex degree below which features for a specific vertex are not generated");
}

void WDK_FeatureGenerator::InitFeatureCache(const GraphClass& aG) {
	mSparseFeatureCache.clear();
	mSparseFeatureFlagCache.clear();
	unsigned vertex_size = aG.VertexSize();
	vector<SVector> init_svector(vertex_size);
	mSparseFeatureCache = init_svector;
	vector<bool> init_bool(vertex_size, false);
	mSparseFeatureFlagCache = init_bool;
}

void WDK_FeatureGenerator::generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList) {
	vector<unsigned> first_endpoint_list = aFirstEndpointList;
	GetFirstEndpoints(aG, first_endpoint_list);
	if (first_endpoint_list.size() == 0)
		throw std::logic_error("ERROR7: Something went wrong: cannot generate features over an empty set of first endpoints!");
	int horizon = max(mDistance, mRadius);
	aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
	if (aG.Check() == false)
		throw logic_error("ERROR12: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

	InitFeatureCache(aG);
	for (unsigned d = 0; d <= mDistance; d++) {
		SVector z;
		for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
			unsigned src_id = first_endpoint_list[i];
			if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
				SVector zv;
				GenerateVertexFeatures(src_id, aG, mRadius, d, zv);
				z.add(zv);
			}
		}
		if (mNormalization)
			z.normalize();
		x.add(z);
	}

	if (mNormalization)
		x.normalize();
	if (mMinKernel)
		ConvertSparseVectorToMinFeatureVector(x);
}

void WDK_FeatureGenerator::generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList) {
	vector<unsigned> first_endpoint_list = aFirstEndpointList;
	GetFirstEndpoints(aG, first_endpoint_list);
	if (first_endpoint_list.size() == 0)
		throw std::logic_error("ERROR8: Something went wrong: cannot generate features over an empty set of first endpoints!");
	int horizon = max(mDistance, mRadius);
	aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
	if (aG.Check() == false)
		throw logic_error("ERROR13: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

	InitFeatureCache(aG);
	for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
		SVector z;
		unsigned src_id = first_endpoint_list[i];
		if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
			for (unsigned d = 0; d <= mDistance; d++) {
				SVector zv;
				GenerateVertexFeatures(src_id, aG, mRadius, d, zv);
				z.add(zv);
			}
		}
		if (mMinKernel)
			ConvertSparseVectorToMinFeatureVector(z);
		x_list.push_back(z);
	}
}

void WDK_FeatureGenerator::GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x) {
	vector<unsigned> endpoint_list(4);
	endpoint_list[0] = aRadius;
	endpoint_list[1] = aDistance;

	unsigned src_code = HashFunc(aG.GetVertexLabelConcatenated(aSrcID));
	vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(aSrcID, aDistance);
	for (unsigned dest_j = 0; dest_j < dest_id_list.size(); dest_j++) {
		unsigned dest_id = dest_id_list[dest_j];
		unsigned dest_code = 0;
		if (aG.GetVertexKernelPoint(dest_id) && aG.GetVertexAlive(dest_id)) { //proceed to extract features only if the *dest* vertex is a kernel point and is alive
			dest_code = HashFunc(aG.GetVertexLabelConcatenated(dest_id));

			//rehash the features in src neighborhood
			endpoint_list[2] = src_code;
			endpoint_list[3] = dest_code;
			unsigned src_code = HashFunc(endpoint_list, mHashBitMask);
			SVector z;
			GenerateVertexFeatures(aSrcID, aG, aRadius, z);
			ReHash(z, src_code);
			x.add(z);

			//rehash the features in dest neighborhood
			endpoint_list[2] = dest_code;
			endpoint_list[3] = src_code;
			unsigned dest_code = HashFunc(endpoint_list, mHashBitMask);
			SVector t;
			GenerateVertexFeatures(dest_j, aG, aRadius, t);
			ReHash(t, dest_code);
			x.add(t);
		}
	}
}

void WDK_FeatureGenerator::ReHash(SVector& x, unsigned aReHashCode) {
	vector<pair<int, double> > vec = x.unpack();
	vector<unsigned> hash_vec(2, 0);
	SVector z;
	for (unsigned i = 0; i < vec.size(); ++i) {
		int key = vec[i].first;
		double val = vec[i].second;
		hash_vec[0] = (unsigned) (key);
		hash_vec[1] = aReHashCode;
		unsigned code = HashFunc(hash_vec, mHashBitMask);
		z.set(code, val);
	}
	x = z;
}

void WDK_FeatureGenerator::GenerateVertexFeatures(unsigned aRootVertexIndex, const GraphClass& aG, unsigned aRadius, SVector& x) {
	if (mSparseFeatureFlagCache[aRootVertexIndex]) {
		x = mSparseFeatureCache[aRootVertexIndex];
	} else {
		unsigned root_degree = aG.VertexAdjacentListSize(aRootVertexIndex);
		if (root_degree <= mLowerVertexDegreeThreshold) { //if the vertex degree is lower than the threshold then just add a dummy feature and return
			SVector z;
			z.set(1, 1);
			x.add(z);
			return;
		}
		//...else
		//extract set of vertices in the ball of radius aRadius
		set<unsigned> ball;
		for (unsigned r = 0; r <= aRadius; r++) {
			vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(aRootVertexIndex, r);
			ball.insert(dest_id_list.begin(), dest_id_list.end());
		}

		//induce the subgraph from the ball and return the new index for the root vertex
		GraphClass gal;
		unsigned root = aG.GetVertexInducedRootedSubGraph(ball, aRootVertexIndex, gal);
		gal.ComputePairwiseDistanceInformation(aRadius * 2);
		string root_label_string = gal.GetVertexLabelConcatenated(root) + ".0";
		unsigned root_label_code = HashFunc(root_label_string);

		vector<unsigned> edge_encoding(4);
		edge_encoding[0] = root_label_code;

		//determine distance of all vertices from root vertex
		vector<unsigned> vertex_distance_list;
		for (unsigned i = 0; i < gal.VertexSize(); ++i) {
			vertex_distance_list.push_back(gal.PairwiseDistance(root, i));
		}
		//for all vertices in ball extract all edges
		for (unsigned u = 0; u < gal.VertexSize(); ++u) {
			string u_label_string = gal.GetVertexLabelConcatenated(u) + "." + stream_cast<string>(vertex_distance_list[u]);
			unsigned u_label_code = HashFunc(u_label_string);
			//get all edges of vertex u
			vector<unsigned> vertex_adjacency_list = gal.GetVertexAdjacentList(u);
			vector<unsigned> edge_adjacency_list = gal.GetEdgeAdjacentList(u);
			for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {
				unsigned v = vertex_adjacency_list[j];
				string v_label_string = gal.GetVertexLabelConcatenated(v) + "." + stream_cast<string>(vertex_distance_list[v]);
				//create a feature as: the root label code + the src vertex label+distance from root code + the dest (farthest) vertex label+distance from root code + the edge label code
				unsigned v_label_code = HashFunc(v_label_string);
				if (u_label_code < v_label_code) {
					edge_encoding[1] = u_label_code;
					edge_encoding[2] = v_label_code;
				} else {
					edge_encoding[1] = v_label_code;
					edge_encoding[2] = u_label_code;
				}
				unsigned e = edge_adjacency_list[j];
				string edge_label = gal.GetEdgeLabelConcatenated(e);
				unsigned hash_edge_label = HashFunc(edge_label);
				edge_encoding[3] = hash_edge_label;
				unsigned hash_edge_encoding = HashFunc(edge_encoding);
				SVector z;
				z.set(hash_edge_encoding, 1);
				x.add(z);
			}
		}
		mSparseFeatureFlagCache[aRootVertexIndex] = true;
		mSparseFeatureCache[aRootVertexIndex] = x;
	}
}

//----------------------------------------------------------------------------------------------------------------------------------------------
ANSPDK_FeatureGenerator::ANSPDK_FeatureGenerator(const std::string& id) :
		NSPDK_FeatureGenerator(id) {
}

void ANSPDK_FeatureGenerator::generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList) {
	vector<unsigned> first_endpoint_list = aFirstEndpointList;
	GetFirstEndpoints(aG, first_endpoint_list);
	if (first_endpoint_list.size() == 0)
		throw std::logic_error("ERROR9: Something went wrong: cannot generate features over an empty set of first endpoints!");
	int horizon = max(mDistance, mRadius);
	aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
	if (aG.Check() == false)
		throw logic_error("ERROR10: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

	InitFeatureCache(aG, mRadius);
	for (unsigned r = 0; r <= mRadius; r++) {
		for (unsigned d = 0; d <= mDistance; d++) {
			SVector z;
			for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
				unsigned src_id = first_endpoint_list[i];
				if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
					GenerateVertexFeatures(src_id, aG, r, d, z);
				} else if (aG.GetVertexAbstraction(src_id) && aG.GetVertexAlive(src_id)) {
					GenerateAbstractVertexFeatures(src_id, aG, r, d, z);
				}
			}
			if (mNormalization)
				z.normalize();
			x.add(z);
		}
	}

	if (mNormalization)
		x.normalize();
	if (mMinKernel)
		ConvertSparseVectorToMinFeatureVector(x);
	if (mDebugVerbosity > 0) {
		cout << x << endl;
		OutputFeatureMap(cout);
		aG.Output(cout);
	}
}

void ANSPDK_FeatureGenerator::generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList) {
	vector<unsigned> first_endpoint_list = aFirstEndpointList;
	GetFirstEndpoints(aG, first_endpoint_list);
	if (first_endpoint_list.size() == 0)
		throw std::logic_error("ERROR6: Something went wrong: cannot generate features over an empty set of first endpoints!");
	int horizon = max(mDistance, mRadius);
	aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
	if (aG.Check() == false)
		throw logic_error("ERROR11: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

	InitFeatureCache(aG, mRadius);
	for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
		SVector z;
		unsigned src_id = first_endpoint_list[i];
		if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
			for (unsigned r = 0; r <= mRadius; r++) {
				for (unsigned d = 0; d <= mDistance; d++) {
					GenerateVertexFeatures(src_id, aG, r, d, z);
				}
			}
		} else if (aG.GetVertexAbstraction(src_id) && aG.GetVertexAlive(src_id)) {
			for (unsigned r = 0; r <= mRadius; r++) {
				for (unsigned d = 0; d <= mDistance; d++) {
					GenerateAbstractVertexFeatures(src_id, aG, r, d, z);
				}
			}
		}
		if (mMinKernel)
			ConvertSparseVectorToMinFeatureVector(z);
		x_list.push_back(z);
	}
	if (mDebugVerbosity > 0) {
		for (unsigned i = 0; i < x_list.size(); ++i)
			cout << i << " " << x_list[i] << endl;
		OutputFeatureMap(cout);
		aG.Output(cout);
	}
}

void ANSPDK_FeatureGenerator::GenerateAbstractVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x) {
	vector<unsigned> endpoint_list(4);
	endpoint_list[0] = aRadius;
	endpoint_list[1] = aDistance;

	//ensure that the src vertex is of the abstraction type
	if (aG.GetVertexAbstraction(aSrcID) == false)
		throw std::logic_error("ERROR2: Something went wrong: expecting an abstraction vertex, but abstraction test fails.");
	//extract all adjacent vertices and partition them into part_of (lower level) and abstraction_of (upper level)
	vector<unsigned> part_of_list;
	vector<unsigned> abstraction_of_list;

	vector<unsigned> vertex_adjacency_list = aG.GetVertexAdjacentList(aSrcID);
	vector<unsigned> edge_adjacency_list = aG.GetEdgeAdjacentList(aSrcID);
	if (vertex_adjacency_list.size() != edge_adjacency_list.size())
		throw std::logic_error("ERROR3: Something went wrong: expecting vertex adjacency list to be the same size as the edge adjacency list.");
	if (vertex_adjacency_list.size() == 0)
		throw std::logic_error("ERROR3b: Something went wrong: expecting a non empty vertex adjacency list for each abstract vertex.");
	for (unsigned i = 0; i < vertex_adjacency_list.size(); ++i) {
		unsigned child_vertex_id = vertex_adjacency_list[i];
		unsigned child_edge_id = edge_adjacency_list[i];
		if (aG.GetEdgePartOf(child_edge_id))
			part_of_list.push_back(child_vertex_id);
		if (aG.GetEdgeAbstractionOf(child_edge_id))
			abstraction_of_list.push_back(child_vertex_id);
	}

	//starting from each upper level vertex find all vertices at distance aDistance
	if (abstraction_of_list.size() == 0)
		throw std::logic_error("ERROR5: Something went wrong: expecting a non empty abstraction_of list.");
	set<unsigned> abstraction_of_set;
	for (unsigned i = 0; i < abstraction_of_list.size(); ++i) {
		unsigned v = abstraction_of_list[i];
		vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(v, aDistance);
		abstraction_of_set.insert(dest_id_list.begin(), dest_id_list.end());
	}
	//starting from each lower level vertex find all vertices at distance aDistance
	if (part_of_list.size() == 0)
		throw std::logic_error("ERROR4: Something went wrong: expecting a non empty part_of list.");
	set<unsigned> part_of_set;
	for (unsigned i = 0; i < part_of_list.size(); ++i) {
		unsigned v = part_of_list[i];
		vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(v, aDistance);
		part_of_set.insert(dest_id_list.begin(), dest_id_list.end());
	}
	//make all possible pairs of distant lower level vertices with distant upper level vertices
	for (set<unsigned>::iterator it = part_of_set.begin(); it != part_of_set.end(); ++it) {
		unsigned part_of_id = *it;
		unsigned part_of_code = GenerateVertexNeighbourhoodHashCode(part_of_id, aG, aRadius);
		endpoint_list[2] = part_of_code;
		for (set<unsigned>::iterator jt = abstraction_of_set.begin(); jt != abstraction_of_set.end(); ++jt) {
			unsigned abstraction_of_id = *jt;
			//build features with one neighborhood graph signature from the lower level and one from the upper
			unsigned abstraction_of_code = GenerateVertexNeighbourhoodHashCode(abstraction_of_id, aG, aRadius);
			endpoint_list[3] = abstraction_of_code;
			unsigned code = HashFunc(endpoint_list, mHashBitMask);
			if (mDebugVerbosity > 0)
				mDebugInfo.StoreFeatureCodeToFeatureInfo(code, endpoint_list);
			SVector z;
			z.set(code, 1);
			x.add(z);
		}
	}
}

//----------------------------------------------------------------------------------------------------------------------------------------------
PBK_FeatureGenerator::PBK_FeatureGenerator(const std::string& id) :
		NSPDK_FeatureGenerator(id), mANSPDK("anspdk" + id), mWDK("wdk" + id) {
}

void PBK_FeatureGenerator::generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList) {
	SVector nspdk_x;
	mANSPDK.generate_feature_vector(aG, nspdk_x, aFirstEndpointList);
	x.add(nspdk_x);
	SVector wdk_x;
	mWDK.generate_feature_vector(aG, wdk_x, aFirstEndpointList);
	x.add(wdk_x);
	if (mNormalization)
		x.normalize();
}

void PBK_FeatureGenerator::generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList) {
	vector<SVector> nspdk_list;
	mANSPDK.generate_vertex_feature_vector(aG, nspdk_list, aFirstEndpointList);
	vector<SVector> wdk_list;
	mWDK.generate_vertex_feature_vector(aG, wdk_list, aFirstEndpointList);
	for (unsigned i = 0; i < nspdk_list.size(); ++i) {
		SVector z;
		z.add(nspdk_list[i]);
		z.add(wdk_list[i]);
		x_list.push_back(z);
	}
}

