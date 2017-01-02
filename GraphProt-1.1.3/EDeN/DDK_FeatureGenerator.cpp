
/*
  * This program is free software; you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation; either version 3 of the License, or
  * (at your option) any later version.
  *
  * Written (W) 2013 Nicolo' Navarin
  */

#include "DDK_FeatureGenerator.h"
 //#include <shogun/features/Labels.h>
#include <climits>
#include <deque>
//#include <boost/container/stable_vector.hpp>
//#include <boost/random/uniform_int.hpp>
//#include <boost/random/random_number_generator.hpp>
//#include <boost/random.hpp>
//#include <boost/generator_iterator.hpp>

//#include <boost/random/variate_generator.hpp>
//#include <boost/generator_iterator.hpp>

// #include <shogun/mathematics/Math.h>

 //TODO eliminare
 double lossy_e=0.001;

 DDkernel_FeatureGenerator::DDkernel_FeatureGenerator()
 : learn_rate(0.1), max_iter(1000), NSPDK_FeatureGenerator("DDK")
 {
 }
/*
 DDkernel_FeatureGenerator::DDkernel_FeatureGenerator(CGraphFeatures* traindat, CLabels* trainlab)
 :  learn_rate(.1), max_iter(1000), NSPDK_FeatureGenerator("DDK")
 {
     //TODO inizializzazione
   //  set_features(traindat);
   //  set_labels(trainlab);
 }*/

 DDkernel_FeatureGenerator::~DDkernel_FeatureGenerator()
 {
	 NSPDK_FeatureGenerator("ddk");
 }
 DDkernel_FeatureGenerator::DDkernel_FeatureGenerator(const std::string& id) :
 		NSPDK_FeatureGenerator(id) {
		new_flag(&mTreeLambda, "mTreeLambda", "(double)\n tree_lambda parameter for DDK kernel");
 }
 void DDkernel_FeatureGenerator::OutputParameters(ostream& out) const {
 	out << "Radius: " << mRadius << endl;
 	out << "Radius2: " << mRadiusTwo << endl;

 	out << "Distance: " << mDistance << endl;
 	out << "Match_Type: " << mMatchType << endl;
 	out << "Hash_Bit_Size: " << mHashBitSize << endl;
 	out << "Hash_Bit_mask: " << mHashBitMask << endl;
 	out << "Min_Kernel: " << mMinKernel << endl;
 	out << "Normalization: " << mNormalization << endl;
 	out << "Vertex Degree Threshold: " << mVertexDegreeThreshold << endl;
 	out << "Debug_Verbosity: " << mDebugVerbosity << endl;
 	out << "Tree_lambda: " << mTreeLambda << endl;

 }

 void DDkernel_FeatureGenerator::generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList) {
 }
	 /*
	 // OutputParameters(cout);
	 vector<GraphClass> graphs;
	 graphs.push_back(aG);
	 CGraphFeatures* data = new CGraphFeatures(graphs);
	// set_lambda(mTreeLambda);
	 //TODO dag_h is now set as the radius parameter r
	// cout<<"dag_h="<<mRadius<<endl;
	 svectorFeatures(data,x,mRadius+1 ,10000000);
	// data->;
	 graphs.clear();
	 delete data;

 }*/



 void DDkernel_FeatureGenerator::generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList) {

 }

 void DDkernel_FeatureGenerator::GenerateAbstractVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x) {

 }

 //---------------------------------------------------------------------------------


//TODO output on svector
 /*
 void DDkernel_FeatureGenerator::svectorFeatures(CGraphFeatures* data,SVector& out_vector, int dag_h, int max_features = 10000000){
      //  vector<mymapForHash*> HashMaps = vector<mymapForHash*>(lhs->get_num_vectors());
     // for di tutti i grafi cosi creo i bigDAG.

 long time = 1.0 / lossy_e;

 vector<double> BER_total;

 OnlineFeatures onlinefeatures(max_features);
      //  OnlineFeatures onlinefeatures_old;

         int64_t numero = 0;



         int training =   data->get_num_vectors();
         string a;
         ofstream output(a.c_str());

         for (int idx_a=0;idx_a<training;idx_a++){

         Graph* a = new Graph() ;
       int32_t size_a = 0;
       ((CGraphFeatures*) data)->get_feature_vector(a,idx_a);

    //   	 try{
       OnlineHashBig2DAGLossy* test = onlinefeatures.generateFeatures(a, dag_h,idx_a,data->get_num_vectors());

       OnlineHashBig2DAGOrdered ordered_map;
    //TODO add features to a OnlineHashBig2DAGOrdered map
       // cout<<"test: "<<test->size();
        BOOST_FOREACH(OnlineHashBig2DAGLossy::value_type local_row, *test) {
           ordered_map[local_row.second.get<2>()]=local_row.second;
       }

    //   cout<<"ordered_map: "<<ordered_map.size();
        //normalize the feature vector


  //      }
        SVector z;

    unsigned int last = 0;
         BOOST_FOREACH(OnlineHashBig2DAGOrdered::value_type i_row,  ordered_map){
           //  cout<<"lambda "<<lambda<<endl;
				z.set((unsigned int)i_row.second.get<2>() % (unsigned int)pow(2,mHashBitSize), pow((double)mTreeLambda, (double)i_row.second.get<0>() / 2.0)* (double)i_row.second.get<3>());

                 ASSERT((unsigned int)i_row.second.get<2>() >= last);
                 last = (unsigned int)i_row.second.get<2>();
         }
      //   if (mNormalization){
        // 				z.normalize();}

		out_vector.add(z);
		//if (mNormalization){
		//		out_vector.normalize();}

     delete a;
            test->clear();
            ordered_map.clear();
         delete test;

     //  	 }catch(shogun::ShogunException e){
     //  		 cout<<e.get_exception_string();

      // 	 }


       }
     //   cout<<endl;
        output.close();


     //  cout<<"FEATURES has been printed"<<endl;



     }





bool GraphSortCriterion (const tuple<double,Graph,long> p1, const tuple<double,Graph,long> p2)
{
    /// a graph is less than another graph


    return p1.get<0>()<p2.get<0>();
}







bool FeaturesSortCriterion (const tuple<double,OnlineHashBig2DAG,long> p1, const tuple<double,OnlineHashBig2DAG,long> p2)
{
    // a graph is less than another graph


    return p1.get<0>()<p2.get<0>();
}

*/

void DDkernel_FeatureGeneratorNew::generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList) {
//	 OnlineFeaturesEDeN of;
//	cout<<"Lambda="<<mTreeLambda<<endl;
//	cout<<"Radius="<<mRadius<<endl;

	 x = generate_feature_vector_core(aG, mRadius);
	 if (mNormalization){
	 				x.normalize();
	 }

}


SVector&  DDkernel_FeatureGeneratorNew::generate_feature_vector_core(const GraphClass& aG , int dag_h){
	//cout<<mTreeLambda<<endl;
	//cout<<"start generation of feature vector.."<<endl;
			SVector* x = new SVector();

		//	cout<<"n of nodes "<<aG.VertexSize()<<endl;
		    for (unsigned i=0; i< aG.VertexSize(); i++){
		    	map<pair<unsigned, unsigned>, int> oSrcDestMaptoDistance;
		    	map<unsigned, vector<unsigned> > MapNodetoParents;
		    	map<unsigned, vector<unsigned> > MapNodetoChildren;

		    	vector<unsigned> TopologicalSort;
		    	int maxLevel = 0;
		    	//aG. // codice calcolo features here
		      	aG.SingleVertexBoundedBreadthFirstVisitTree(i,dag_h, oSrcDestMaptoDistance,  MapNodetoParents,MapNodetoChildren, TopologicalSort, maxLevel);
		//    cout<<"MaxLevel: "<<maxLevel<<endl;
		      	// 	cout<<"Breadth first visit ended.."<<endl;
		      	// I have the DAG for the vertex i. Need topological sort.
		      	//map every node to him productions
		    	map<unsigned, vector<unsigned> > MapNodeToProductionsID;
		    	map<unsigned, vector<int> > MapNodetoFrequencies;

		    	map<unsigned, int > MapProductionIDtoSize;

		  //  	cout<<"topological sort size "<<TopologicalSort.size()<<endl;

		      	//maybe I meed the maximum depth (distance)
		      	for (int ii=TopologicalSort.size()-1; ii>=0; ii--){
		      		// cout<<ii<<endl;
		      		//calculate child size
			    	vector<unsigned> vertex_adjacency_list = MapNodetoChildren[TopologicalSort[ii]];
		      		int max_child_heigth=0;
		        	for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {
		    			int child_vertex_id= vertex_adjacency_list[j];
		    			int child_heigth=MapNodeToProductionsID[child_vertex_id].size();
		    			if(child_heigth > max_child_heigth){
		    				max_child_heigth = child_heigth;
		    			}
		        	}

		      		for (int depth=0; depth<= max_child_heigth; ++depth){ //era maxLevel
		      			unsigned hash_subgraph_code = 1;
		      		//depth 0; only label
		      			if(depth==0){
		      		//		cout<<"depth=0"<<endl;
		      				unsigned enc= Radius0RootedGraphCanonicalFormEncoding(TopologicalSort[ii], aG);
		      		MapNodeToProductionsID[TopologicalSort[ii]].push_back(enc);
		      		// scorro i figli e creo etichetta tramite funzione di hash
		      		// depth > 0; the length of the children ID list
		      		//imposta valore feature a 1 (TODO lambda)
		      		//							level(o-mRadius)
			    	int frequency = 0;
			    	if(max_child_heigth==0){
		      		//leaf node
			    	frequency=  maxLevel - oSrcDestMaptoDistance[make_pair(i,TopologicalSort[ii])];
			    	}
			    	SVector z;
		      		z.set(enc, (double)(frequency+1.0) * sqrt(mTreeLambda));
			    //	cout<<"Feature: "<<aG.GetVertexLabelConcatenated(TopologicalSort[ii])<<" freq:"<<(double)(frequency+1.0) * sqrt(mTreeLambda)<<" code: "<<enc<<endl;
			    	x->add(z);
			    	MapNodetoFrequencies[TopologicalSort[ii]].push_back(frequency);
		      		MapProductionIDtoSize[enc]= 1; //sqrt(mTreeLambda);
		      			}
		      			else{
		      	//			cout<<"depth="<<depth<<endl;

		      			int size = 0;
		      		string encoding;
		      		encoding = aG.GetVertexLabelConcatenated(TopologicalSort[ii]);
			    	vector<unsigned> vertex_adjacency_list = MapNodetoChildren[TopologicalSort[ii]];
		    		vector<unsigned>  vertex_label_id_list;
		    		int min_freq_children =INT_MAX;
			    	for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {

			    		unsigned child_vertex_id = vertex_adjacency_list[j];
			    		int size_map =  MapNodeToProductionsID[child_vertex_id].size();//min (depth-1) e dim
			    		unsigned child_hash = MapNodeToProductionsID[child_vertex_id][min(size_map,depth)-1];// id hash del figlio;
			    		int freq_child = MapNodetoFrequencies[child_vertex_id][min(size_map,depth)-1];
			    		if (freq_child < min_freq_children){
			    			min_freq_children = freq_child;
			    		}
			    		vertex_label_id_list.push_back(child_hash);
			    		size += MapProductionIDtoSize[child_hash];
			    	}

			    	// vertex_label_id_list contiene gli hash delle prod dei figli. Va ordinato e usato x generare l hash attuale
			    	sort(vertex_label_id_list.begin(), vertex_label_id_list.end());

			    	if (vertex_label_id_list.size() > 0){
			    		std::stringstream ss;
			    		ss <<":"<<vertex_label_id_list[0];
			    		encoding += ss.str() ;
			    	}//originale senza .
			    	for (unsigned i = 1; i < vertex_label_id_list.size(); i++){
			    		std::stringstream ss;
			    		ss << vertex_label_id_list[i];
			    		encoding += "." + ss.str() ;
			    	}
			    //DEBUG	cout<<"encoding "<<encoding<<endl;
			    	//calculate frequency

			    	hash_subgraph_code = HashFunc(encoding);
			    	MapNodeToProductionsID[TopologicalSort[ii]].push_back(hash_subgraph_code);
			    	size += 1; // current node
			    	MapProductionIDtoSize[hash_subgraph_code]= size;
			    	int frequency = min_freq_children;
			    	MapNodetoFrequencies[TopologicalSort[ii]].push_back(frequency);

			    	SVector z;
			   // 	cout<<"Feature: "<<encoding<<" freq:"<<(double)(frequency+1.0) * sqrt(pow(mTreeLambda,size))<<" code: "<<hash_subgraph_code<<endl;
			    	z.set(hash_subgraph_code,(double)(frequency+1.0) * sqrt(pow(mTreeLambda,size)));
			    	x->add(z);

			    //	cout<<"weight "<<sqrt(pow(mTreeLambda,size))<<endl;

		      			}

		      		}
		      	}


		    }
		    return *x;

	 }

DDkernel_FeatureGeneratorNew::DDkernel_FeatureGeneratorNew(const std::string& id) :
		DDkernel_FeatureGenerator(id) {
}

NSDDkernel_FeatureGenerator::NSDDkernel_FeatureGenerator(const std::string& id) :
		DDkernel_FeatureGeneratorNew(id) {
	mRadiusTwo = 0;
			new_flag(&mRadiusTwo, "mRadiusTwo", "(double)\n second radius parameter for NSDDK kernel");

}

void NSDDkernel_FeatureGenerator::generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList){
//cout<<mRadiusTwo<<endl;
	 x = generate_feature_vector_core(aG, mRadius);
	//x = generate_feature_vector_core_DISTANCE(aG, mRadius); //TEST distanze

	// if (mNormalization){
	 //				x.normalize();
	// }
}


SVector&   NSDDkernel_FeatureGenerator::generate_feature_vector_core(const GraphClass& aG , int dag_h){
   // test rNSPDK 1/2 d=1/2
	GraphClass ag2 = aG;
	SVector* out = new SVector();
//	cout<<"Generate feature vector core"<<endl;
	for (unsigned i=0; i< aG.VertexSize(); i++){
	//	cout<<"i:"<<i<<endl;

		//SVector DD= generate_vertex_feature_vector(i,ag2 ,dag_h);
		vector<SVector> features; //(dag_h+1);
		SVector DD= generate_vertex_feature_vector_DIVIDED(i,ag2 ,mRadius,features);
	//	cout<<"got DD vector"<<endl;
		//cout<<features.size();
	//	out->add(DD);
		//NSPDK
	//	cout<<"start generation of NSPKD feature vector.. mDistance"<<mDistance<<endl;

		//NSPDK_FeatureGenerator::GenerateVertexFeatures(i,aG,NSPDK);
	//	NSPDK_FeatureGenerator::generate_feature_vector(aG,NSPDK);
		vector<unsigned> first_endpoint_list= vector<unsigned>();
		first_endpoint_list.push_back(i);
			//GetFirstEndpoints(aG, first_endpoint_list);

		int horizon=max(mDistance, mRadiusTwo); //   max(mDistance, mRadius);
			aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
			if (aG.Check() == false)
				throw logic_error("ERROR10: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

			InitFeatureCache(aG, mRadiusTwo ); //mRadius
			for (unsigned r = 0; r <= mRadius ; r++) { //max(1.0,mRadius / 2.0)
				for (unsigned r2 = 0; r2 <= mRadiusTwo ; r2++) {
				for (unsigned d = 0; d <= mDistance; d++) { // era 0
					SVector z;
					for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
						unsigned src_id = first_endpoint_list[i];
						if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
							SVector zv;

							vector<unsigned> endpoint_list(4);
							endpoint_list[0] = r2;
							endpoint_list[1] = d;

						//	unsigned src_code = GenerateVertexNeighbourhoodHashCode(aSrcID, aG, aRadius);
				//	TODO	for(unsigned src_code in DD)

							vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(i, d);
							for (unsigned dest_j = 0; dest_j < dest_id_list.size(); dest_j++) {
								unsigned dest_id = dest_id_list[dest_j];
								unsigned dest_code = 0;
								if (aG.GetVertexKernelPoint(dest_id) && aG.GetVertexAlive(dest_id)) { //proceed to extract features only if the *dest* vertex is a kernel point and is alive
									dest_code = GenerateVertexNeighbourhoodHashCode(dest_id, aG, r2); //
									// first is the DDK feature, second NSPDK feature
						//  for(const SVector::Pair *p = DD; p->i>=0; p++) {
						 for(const SVector::Pair *p = features[r]; p->i>=0; p++) {

									endpoint_list[2] = p->i;
										//TODO scorrere le feature di DD
										endpoint_list[3] = dest_code;
									unsigned code = HashFunc(endpoint_list, mHashBitMask);

									SVector z;
									z.set(code, p->v);
									zv.add(z);
								}
								}
							}



							z.add(zv);
						}
						}

					if (mNormalization)
						z.normalize();
				out->add(z);
				} //d
			} //r
			} //r2







	//	cout<<"features size DD:"<<DD.sparse_size()<<" out:"<<out->sparse_size()<<endl;


}
	return *out;
}

SVector&   NSDDkernel_FeatureGenerator::generate_feature_vector_core_DISTANCE(const GraphClass& aG , int dag_h){
   // test rNSPDK 1/2 d=1/2
	GraphClass ag2 = aG;
	SVector* out = new SVector();
//	cout<<"Generate feature vector core"<<endl;
	for (unsigned i=0; i< aG.VertexSize(); i++){
	//	cout<<"i:"<<i<<endl;

		//SVector DD= generate_vertex_feature_vector(i,ag2 ,dag_h);
		vector< vector< SVector> > features; //(dag_h+1);
		SVector DD= generate_vertex_feature_vector_DIVIDED_DISTANCE(i,ag2 ,dag_h,features);
		//cout<<features.size();
	//	out->add(DD);
		//NSPDK
	//	cout<<"start generation of NSPKD feature vector.. mDistance"<<mDistance<<endl;

		//NSPDK_FeatureGenerator::GenerateVertexFeatures(i,aG,NSPDK);
	//	NSPDK_FeatureGenerator::generate_feature_vector(aG,NSPDK);
		vector<unsigned> first_endpoint_list= vector<unsigned>();
		first_endpoint_list.push_back(i);
			//GetFirstEndpoints(aG, first_endpoint_list);

		int horizon=max(mDistance, mRadiusTwo); //   max(mDistance, mRadius);
			aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
			if (aG.Check() == false)
				throw logic_error("ERROR10: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

			InitFeatureCache(aG, mRadiusTwo ); //mRadius
			for (unsigned r = 0; r <= mRadius ; r++) { //max(1.0,mRadius / 2.0)
				for (unsigned r2 = 0; r2 <= mRadiusTwo ; r2++) {
				for (unsigned d = 0; d <= mDistance; d++) { // era 0
					SVector z;
					for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
						unsigned src_id = first_endpoint_list[i];
						if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
							SVector zv;

							vector<unsigned> endpoint_list(5);
							endpoint_list[0] = r2;
							endpoint_list[1] = d;

						//	unsigned src_code = GenerateVertexNeighbourhoodHashCode(aSrcID, aG, aRadius);
				//	TODO	for(unsigned src_code in DD)

							vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(i, d);
							for (unsigned dest_j = 0; dest_j < dest_id_list.size(); dest_j++) {
								unsigned dest_id = dest_id_list[dest_j];
								unsigned dest_code = 0;
								if (aG.GetVertexKernelPoint(dest_id) && aG.GetVertexAlive(dest_id)) { //proceed to extract features only if the *dest* vertex is a kernel point and is alive
									dest_code = GenerateVertexNeighbourhoodHashCode(dest_id, aG, r2); //
									// first is the DDK feature, second NSPDK feature
						//  for(const SVector::Pair *p = DD; p->i>=0; p++) {
									for(int distance=0; distance<features[r].size(); distance++){ //dag_h+1
										endpoint_list[4]=distance; //from root in DD
						 for(const SVector::Pair *p = features[r][distance]; p->i>=0; p++) {

									endpoint_list[2] = p->i;
										//TODO scorrere le feature di DD
										endpoint_list[3] = dest_code;
										//TEST output
								//		cout<<"R "<<endpoint_list[0]<<" d "<<endpoint_list[1]<<" DD "<<endpoint_list[2]<<" NSPDK "<<endpoint_list[3]<<" dRoot "<<endpoint_list[4]<<endl;
									unsigned code = HashFunc(endpoint_list, mHashBitMask);

									SVector z;
									z.set(code, p->v);
									zv.add(z);
								}
								}
								}
							}



							z.add(zv);
						}
						}

					if (mNormalization)
						z.normalize();
				out->add(z);
				} //d
			} //r
			}//r2







	//	cout<<"features size DD:"<<DD.sparse_size()<<" out:"<<out->sparse_size()<<endl;


}
	return *out;
}

SVector&   NSDDkernel_FeatureGenerator::generate_vertex_feature_vector(unsigned i, const GraphClass& aG , int dag_h){
	//cout<<"start generation of feature vector.."<<endl;
			SVector* x = new SVector();

		//	cout<<"n of nodes "<<aG.VertexSize()<<endl;
		    assert(i< aG.VertexSize());
		    	map<pair<unsigned, unsigned>, int> oSrcDestMaptoDistance;
		    	map<unsigned, vector<unsigned> > MapNodetoParents;
		    	map<unsigned, vector<unsigned> > MapNodetoChildren;

		    	vector<unsigned> TopologicalSort;
		    	int maxLevel = 0;
		    	//aG. // codice calcolo features here
		      	aG.SingleVertexBoundedBreadthFirstVisitTree(i,dag_h, oSrcDestMaptoDistance,  MapNodetoParents,MapNodetoChildren, TopologicalSort, maxLevel);
		//    cout<<"MaxLevel: "<<maxLevel<<endl;
		      	// 	cout<<"Breadth first visit ended.."<<endl;
		      	// I have the DAG for the vertex i. Need topological sort.
		      	//map every node to him productions
		    	map<unsigned, vector<unsigned> > MapNodeToProductionsID;
		    	map<unsigned, vector<int> > MapNodetoFrequencies;

		    	map<unsigned, int > MapProductionIDtoSize;

		  //  	cout<<"topological sort size "<<TopologicalSort.size()<<endl;

		      	//maybe I meed the maximum depth (distance)
		      	for (int ii=TopologicalSort.size()-1; ii>=0; ii--){
		      		// cout<<ii<<endl;
		      		//calculate child size
			    	vector<unsigned> vertex_adjacency_list = MapNodetoChildren[TopologicalSort[ii]];
		      		int max_child_heigth=0;
		        	for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {
		    			int child_vertex_id= vertex_adjacency_list[j];
		    			int child_heigth=MapNodeToProductionsID[child_vertex_id].size();
		    			if(child_heigth > max_child_heigth){
		    				max_child_heigth = child_heigth;
		    			}
		        	}

		      		for (int depth=0; depth<= max_child_heigth; ++depth){ //era maxLevel
		      			unsigned hash_subgraph_code = 1;
		      		//depth 0; only label
		      			if(depth==0){
		      		//		cout<<"depth=0"<<endl;
		      				unsigned enc= Radius0RootedGraphCanonicalFormEncoding(TopologicalSort[ii], aG);
		      		MapNodeToProductionsID[TopologicalSort[ii]].push_back(enc);
		      		// scorro i figli e creo etichetta tramite funzione di hash
		      		// depth > 0; the length of the children ID list
		      		//imposta valore feature a 1 (TODO lambda)
		      		//							level(o-mRadius)
			    	int frequency = 0;
			    	if(max_child_heigth==0){
		      		//leaf node
			    	frequency=  maxLevel - oSrcDestMaptoDistance[make_pair(i,TopologicalSort[ii])];
			    	}
			    	SVector z;
		      		z.set(enc, (double)(frequency+1.0) * sqrt(mTreeLambda));
			    //	cout<<"Feature: "<<aG.GetVertexLabelConcatenated(TopologicalSort[ii])<<" freq:"<<(double)(frequency+1.0) * sqrt(mTreeLambda)<<" code: "<<enc<<endl;
			    	x->add(z);
			    	MapNodetoFrequencies[TopologicalSort[ii]].push_back(frequency);
		      		MapProductionIDtoSize[enc]= 1; //sqrt(mTreeLambda);
		      			}
		      			else{
		      	//			cout<<"depth="<<depth<<endl;

		      			int size = 0;
		      		string encoding;
		      		encoding = aG.GetVertexLabelConcatenated(TopologicalSort[ii]);
			    	vector<unsigned> vertex_adjacency_list = MapNodetoChildren[TopologicalSort[ii]];
		    		vector<unsigned>  vertex_label_id_list;
		    		int min_freq_children =INT_MAX;
			    	for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {

			    		unsigned child_vertex_id = vertex_adjacency_list[j];
			    		int size_map =  MapNodeToProductionsID[child_vertex_id].size();//min (depth-1) e dim
			    		unsigned child_hash = MapNodeToProductionsID[child_vertex_id][min(size_map,depth)-1];// id hash del figlio;
			    		int freq_child = MapNodetoFrequencies[child_vertex_id][min(size_map,depth)-1];
			    		if (freq_child < min_freq_children){
			    			min_freq_children = freq_child;
			    		}
			    		vertex_label_id_list.push_back(child_hash);
			    		size += MapProductionIDtoSize[child_hash];
			    	}

			    	// vertex_label_id_list contiene gli hash delle prod dei figli. Va ordinato e usato x generare l hash attuale
			    	sort(vertex_label_id_list.begin(), vertex_label_id_list.end());

			    	if (vertex_label_id_list.size() > 0){
			    		std::stringstream ss;
			    		ss <<":"<<vertex_label_id_list[0];
			    		encoding += ss.str() ;
			    	}//originale senza .
			    	for (unsigned i = 1; i < vertex_label_id_list.size(); i++){
			    		std::stringstream ss;
			    		ss << vertex_label_id_list[i];
			    		encoding += "." + ss.str() ;
			    	}
			    //DEBUG	cout<<"encoding "<<encoding<<endl;
			    	//calculate frequency

			    	hash_subgraph_code = HashFunc(encoding);
			    	MapNodeToProductionsID[TopologicalSort[ii]].push_back(hash_subgraph_code);
			    	size += 1; // current node
			    	MapProductionIDtoSize[hash_subgraph_code]= size;
			    	int frequency = min_freq_children;
			    	MapNodetoFrequencies[TopologicalSort[ii]].push_back(frequency);

			    	SVector z;
			   // 	cout<<"Feature: "<<encoding<<" freq:"<<(double)(frequency+1.0) * sqrt(pow(mTreeLambda,size))<<" code: "<<hash_subgraph_code<<endl;
			    	z.set(hash_subgraph_code,(double)(frequency+1.0) * sqrt(pow(mTreeLambda,size)));
			    	x->add(z);

			    //	cout<<"weight "<<sqrt(pow(mTreeLambda,size))<<endl;

		      			}

		      		}
		      	}



		    return *x;

	 }


SVector&   NSDDkernel_FeatureGenerator::generate_vertex_feature_vector_DIVIDED(unsigned i, const GraphClass& aG , int dag_h, vector<SVector> & features_grouped_by_depth){
	//cout<<"start generation of feature vector.."<<endl;
			SVector* x = new SVector();
			features_grouped_by_depth= vector<SVector>(dag_h+1);
		/*	for(int i=0; i<dag_h;i++){
				features_grouped_by_depth.push_back(SVector(0));

			}*/
		//	cout<<"n of nodes "<<aG.VertexSize()<<endl;
		    assert(i< aG.VertexSize());
		    	map<pair<unsigned, unsigned>, int> oSrcDestMaptoDistance;
		    	map<unsigned, vector<unsigned> > MapNodetoParents;
		    	map<unsigned, vector<unsigned> > MapNodetoChildren;

		    	vector<unsigned> TopologicalSort;
		    	int maxLevel = 0;
		    	//aG. // codice calcolo features here
		      	aG.SingleVertexBoundedBreadthFirstVisitTree(i,dag_h, oSrcDestMaptoDistance,  MapNodetoParents,MapNodetoChildren, TopologicalSort, maxLevel);
		//    cout<<"MaxLevel: "<<maxLevel<<endl;
		      	// 	cout<<"Breadth first visit ended.."<<endl;
		      	// I have the DAG for the vertex i. Need topological sort.
		      	//map every node to him productions
		    	map<unsigned, vector<unsigned> > MapNodeToProductionsID;
		    	map<unsigned, vector<int> > MapNodetoFrequencies;

		    	map<unsigned, int > MapProductionIDtoSize;

		  //  	cout<<"topological sort size "<<TopologicalSort.size()<<endl;

		      	//maybe I meed the maximum depth (distance)
		      	for (int ii=TopologicalSort.size()-1; ii>=0; ii--){
		      		// cout<<ii<<endl;
		      		//calculate child size
			    	vector<unsigned> vertex_adjacency_list = MapNodetoChildren[TopologicalSort[ii]];
		      		int max_child_heigth=0;
		        	for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {
		    			int child_vertex_id= vertex_adjacency_list[j];
		    			int child_heigth=MapNodeToProductionsID[child_vertex_id].size();
		    			if(child_heigth > max_child_heigth){
		    				max_child_heigth = child_heigth;
		    			}
		        	}

		      		for (int depth=0; depth<= max_child_heigth; ++depth){ //era maxLevel
		      			unsigned hash_subgraph_code = 1;
		      		//depth 0; only label
		      			if(depth==0){
		      		//		cout<<"depth=0"<<endl;
		      				unsigned enc= Radius0RootedGraphCanonicalFormEncoding(TopologicalSort[ii], aG);
		      		MapNodeToProductionsID[TopologicalSort[ii]].push_back(enc);
		      		// scorro i figli e creo etichetta tramite funzione di hash
		      		// depth > 0; the length of the children ID list
		      		//imposta valore feature a 1 (TODO lambda)
		      		//							level(o-mRadius)
			    	int frequency = 0;
			    	if(max_child_heigth==0){
		      		//leaf node
			    	frequency=  maxLevel - oSrcDestMaptoDistance[make_pair(i,TopologicalSort[ii])];
			    	}
			    	SVector z;
		      		z.set(enc, (double)(frequency+1.0) * sqrt(mTreeLambda));
			    //	cout<<"Feature: "<<aG.GetVertexLabelConcatenated(TopologicalSort[ii])<<" freq:"<<(double)(frequency+1.0) * sqrt(mTreeLambda)<<" code: "<<enc<<endl;
			    	x->add(z);
			    	features_grouped_by_depth[depth].add(z);
			    	MapNodetoFrequencies[TopologicalSort[ii]].push_back(frequency);
		      		MapProductionIDtoSize[enc]= 1; //sqrt(mTreeLambda);
		      			}
		      			else{
		      	//			cout<<"depth="<<depth<<endl;

		      			int size = 0;
		      		string encoding;
		      		encoding = aG.GetVertexLabelConcatenated(TopologicalSort[ii]);
			    	vector<unsigned> vertex_adjacency_list = MapNodetoChildren[TopologicalSort[ii]];
		    		vector<unsigned>  vertex_label_id_list;
		    		int min_freq_children =INT_MAX;
			    	for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {

			    		unsigned child_vertex_id = vertex_adjacency_list[j];
			    		int size_map =  MapNodeToProductionsID[child_vertex_id].size();//min (depth-1) e dim
			    		unsigned child_hash = MapNodeToProductionsID[child_vertex_id][min(size_map,depth)-1];// id hash del figlio;
			    		int freq_child = MapNodetoFrequencies[child_vertex_id][min(size_map,depth)-1];
			    		if (freq_child < min_freq_children){
			    			min_freq_children = freq_child;
			    		}
			    		vertex_label_id_list.push_back(child_hash);
			    		size += MapProductionIDtoSize[child_hash];
			    	}

			    	// vertex_label_id_list contiene gli hash delle prod dei figli. Va ordinato e usato x generare l hash attuale
			    	sort(vertex_label_id_list.begin(), vertex_label_id_list.end());

			    	if (vertex_label_id_list.size() > 0){
			    		std::stringstream ss;
			    		ss <<":"<<vertex_label_id_list[0];
			    		encoding += ss.str() ;
			    	}//originale senza .
			    	for (unsigned i = 1; i < vertex_label_id_list.size(); i++){
			    		std::stringstream ss;
			    		ss << vertex_label_id_list[i];
			    		encoding += "." + ss.str() ;
			    	}
			    //DEBUG	cout<<"encoding "<<encoding<<endl;
			    	//calculate frequency

			    	hash_subgraph_code = HashFunc(encoding);
			    	MapNodeToProductionsID[TopologicalSort[ii]].push_back(hash_subgraph_code);
			    	size += 1; // current node
			    	MapProductionIDtoSize[hash_subgraph_code]= size;
			    	int frequency = min_freq_children;
			    	MapNodetoFrequencies[TopologicalSort[ii]].push_back(frequency);

			    	SVector z;
			   // 	cout<<"Feature: "<<encoding<<" freq:"<<(double)(frequency+1.0) * sqrt(pow(mTreeLambda,size))<<" code: "<<hash_subgraph_code<<endl;
			    	z.set(hash_subgraph_code,(double)(frequency+1.0) * sqrt(pow(mTreeLambda,size)));
			    	x->add(z);
			    	features_grouped_by_depth[depth].add(z);
			    //	cout<<"weight "<<sqrt(pow(mTreeLambda,size))<<endl;

		      			}

		      		}
		      	}



		    return *x;

	 }


SVector&   NSDDkernel_FeatureGenerator::generate_vertex_feature_vector_DIVIDED_DISTANCE(unsigned i, const GraphClass& aG , int dag_h, vector<vector<SVector> > & features_grouped_by_depth){
	//cout<<"start generation of feature vector.."<<endl;
			SVector* x = new SVector();
			features_grouped_by_depth= vector<vector<SVector> >(dag_h+1);


		//	cout<<"n of nodes "<<aG.VertexSize()<<endl;
		    assert(i< aG.VertexSize());
		    	map<pair<unsigned, unsigned>, int> oSrcDestMaptoDistance;
		    	map<unsigned, vector<unsigned> > MapNodetoParents;
		    	map<unsigned, vector<unsigned> > MapNodetoChildren;

		    	vector<unsigned> TopologicalSort;
		    	int maxLevel = 0;
		    	//aG. // codice calcolo features here
		      	aG.SingleVertexBoundedBreadthFirstVisitTree(i,dag_h, oSrcDestMaptoDistance,  MapNodetoParents,MapNodetoChildren, TopologicalSort, maxLevel);
		//    cout<<"MaxLevel: "<<maxLevel<<endl;
				for(int t=0; t<maxLevel+1;t++){
					features_grouped_by_depth[t] =(vector<SVector>(dag_h+1));

				}

		      	// 	cout<<"Breadth first visit ended.."<<endl;
		      	// I have the DAG for the vertex i. Need topological sort.
		      	//map every node to him productions
		    	map<unsigned, vector<unsigned> > MapNodeToProductionsID;
		    	map<unsigned, vector<int> > MapNodetoFrequencies;

		    	map<unsigned, int > MapProductionIDtoSize;

		  //  	cout<<"topological sort size "<<TopologicalSort.size()<<endl;

		      	//maybe I meed the maximum depth (distance)
		      	for (int ii=TopologicalSort.size()-1; ii>=0; ii--){
		      		// cout<<ii<<endl;
		      		//calculate child size
			    	vector<unsigned> vertex_adjacency_list = MapNodetoChildren[TopologicalSort[ii]];
		      		int max_child_heigth=0;
		        	for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {
		    			int child_vertex_id= vertex_adjacency_list[j];
		    			int child_heigth=MapNodeToProductionsID[child_vertex_id].size();
		    			if(child_heigth > max_child_heigth){
		    				max_child_heigth = child_heigth;
		    			}
		        	}

		      		for (int depth=0; depth<= max_child_heigth; ++depth){ //era maxLevel
		      			unsigned hash_subgraph_code = 1;
		      		//depth 0; only label
		      			if(depth==0){
		      		//		cout<<"depth=0"<<endl;
		      				unsigned enc= Radius0RootedGraphCanonicalFormEncoding(TopologicalSort[ii], aG);
		      		MapNodeToProductionsID[TopologicalSort[ii]].push_back(enc);
		      		// scorro i figli e creo etichetta tramite funzione di hash
		      		// depth > 0; the length of the children ID list
		      		//imposta valore feature a 1 (TODO lambda)
		      		//							level(o-mRadius)
			    	int frequency = 0;
			    	if(max_child_heigth==0){
		      		//leaf node
			    	frequency=  maxLevel - oSrcDestMaptoDistance[make_pair(i,TopologicalSort[ii])];
			    	}
			    	SVector z;
		      		z.set(enc, (double)(frequency+1.0) * sqrt(mTreeLambda));
			    //	cout<<"Feature: "<<aG.GetVertexLabelConcatenated(TopologicalSort[ii])<<" freq:"<<(double)(frequency+1.0) * sqrt(mTreeLambda)<<" code: "<<enc<<endl;
			    	x->add(z);
			    	features_grouped_by_depth[depth][oSrcDestMaptoDistance[make_pair(i,TopologicalSort[ii])]].add(z);
			    	MapNodetoFrequencies[TopologicalSort[ii]].push_back(frequency);
		      		MapProductionIDtoSize[enc]= 1; //sqrt(mTreeLambda);
		      			}
		      			else{
		      	//			cout<<"depth="<<depth<<endl;

		      			int size = 0;
		      		string encoding;
		      		encoding = aG.GetVertexLabelConcatenated(TopologicalSort[ii]);
			    	vector<unsigned> vertex_adjacency_list = MapNodetoChildren[TopologicalSort[ii]];
		    		vector<unsigned>  vertex_label_id_list;
		    		int min_freq_children =INT_MAX;
			    	for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {

			    		unsigned child_vertex_id = vertex_adjacency_list[j];
			    		int size_map =  MapNodeToProductionsID[child_vertex_id].size();//min (depth-1) e dim
			    		unsigned child_hash = MapNodeToProductionsID[child_vertex_id][min(size_map,depth)-1];// id hash del figlio;
			    		int freq_child = MapNodetoFrequencies[child_vertex_id][min(size_map,depth)-1];
			    		if (freq_child < min_freq_children){
			    			min_freq_children = freq_child;
			    		}
			    		vertex_label_id_list.push_back(child_hash);
			    		size += MapProductionIDtoSize[child_hash];
			    	}

			    	// vertex_label_id_list contiene gli hash delle prod dei figli. Va ordinato e usato x generare l hash attuale
			    	sort(vertex_label_id_list.begin(), vertex_label_id_list.end());

			    	if (vertex_label_id_list.size() > 0){
			    		std::stringstream ss;
			    		ss <<":"<<vertex_label_id_list[0];
			    		encoding += ss.str() ;
			    	}//originale senza .
			    	for (unsigned id = 1; id < vertex_label_id_list.size(); id++){
			    		std::stringstream ss;
			    		ss << vertex_label_id_list[id];
			    		encoding += "." + ss.str() ;
			    	}
			    //DEBUG	cout<<"encoding "<<encoding<<endl;
			    	//calculate frequency

			    	hash_subgraph_code = HashFunc(encoding);
			    	MapNodeToProductionsID[TopologicalSort[ii]].push_back(hash_subgraph_code);
			    	size += 1; // current node
			    	MapProductionIDtoSize[hash_subgraph_code]= size;
			    	int frequency = min_freq_children;
			    	MapNodetoFrequencies[TopologicalSort[ii]].push_back(frequency);

			    	SVector z;
			   // 	cout<<"Feature: "<<encoding<<" freq:"<<(double)(frequency+1.0) * sqrt(pow(mTreeLambda,size))<<" code: "<<hash_subgraph_code<<endl;
			    	z.set(hash_subgraph_code,(double)(frequency+1.0) * sqrt(pow(mTreeLambda,size)));
			    	x->add(z);
			    	features_grouped_by_depth[depth][oSrcDestMaptoDistance[make_pair(i,TopologicalSort[ii])]].add(z); //TODO distance from root
			    //	cout<<"weight "<<sqrt(pow(mTreeLambda,size))<<endl;

		      			}

		      		}
		      	}



		    return *x;

	 }

//Kernel su dati con astrazione (nodi di relazione)
ANSDDkernel_FeatureGenerator::ANSDDkernel_FeatureGenerator(const std::string& id) :
		NSDDkernel_FeatureGenerator(id) {
		//new_flag(&mTreeLambda, "mTreeLambda", "(double)\n tree_lambda parameter for DDK kernel");
}


void ANSDDkernel_FeatureGenerator::generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList) {
	//cout<<"generate_vertex_feature_vector.... OK for abstract"<<endl;
	vector<unsigned> first_endpoint_list = aFirstEndpointList;
	GetFirstEndpoints(aG, first_endpoint_list);
	if (first_endpoint_list.size() == 0)
		throw std::logic_error("ERROR6: Something went wrong: cannot generate features over an empty set of first endpoints!");
	int horizon = max(mDistance, mRadius);
	aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
	if (aG.Check() == false)
		throw logic_error("ERROR11: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

	InitFeatureCache(aG, mRadiusTwo);
	for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
		SVector z;
		unsigned src_id = first_endpoint_list[i];
		//InitFeatureCache(aG, mRadiusTwo);

		if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
			//for (unsigned r = 0; r <= mRadius; r++) {
			//	for (unsigned d = 0; d <= mDistance; d++) {
				//	GenerateVertexFeatures(src_id, aG, r, d, z);
			GenerateVertexFeatures(src_id, aG, mRadius, mDistance, z);
				//}
			//}
		} else if (aG.GetVertexAbstraction(src_id) && aG.GetVertexAlive(src_id)) {
			//for (unsigned r = 0; r <= mRadius; r++) {
			//	for (unsigned d = 0; d <= mDistance; d++) {
				//	GenerateAbstractVertexFeatures(src_id, aG, r, d, z);
			GenerateAbstractVertexFeatures(src_id, aG, mRadius, mDistance, z);

			//	}
			//}
		}
		if (mMinKernel)
			ConvertSparseVectorToMinFeatureVector(z);
		//x_list.push_back(z);
		x.add(z);
	}
	if (mDebugVerbosity > 0) {

		OutputFeatureMap(cout);
		aG.Output(cout);
	}
}
void ANSDDkernel_FeatureGenerator::GenerateVertexFeatures(unsigned aRootVertexIndex, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x) {
	x= generateVertexFeatures_DISTANCE(aRootVertexIndex, aG ,aRadius);
}

SVector&   ANSDDkernel_FeatureGenerator::generateVertexFeatures_DISTANCE(unsigned i, const GraphClass& aG , int dag_h){
   // test rNSPDK 1/2 d=1/2
	GraphClass ag2 = aG;
	SVector* out = new SVector();
//	cout<<"Generate feature vector core"<<endl;
	//	cout<<"i:"<<i<<endl;

		//SVector DD= generate_vertex_feature_vector(i,ag2 ,dag_h);
		vector< vector< SVector> > features; //(dag_h+1);
		SVector DD= generate_vertex_feature_vector_DIVIDED_DISTANCE(i,ag2 ,dag_h,features);
		//cout<<features.size();
	//	out->add(DD);
		//NSPDK
	//	cout<<"start generation of NSPKD feature vector.. mDistance"<<mDistance<<endl;

		//NSPDK_FeatureGenerator::GenerateVertexFeatures(i,aG,NSPDK);
	//	NSPDK_FeatureGenerator::generate_feature_vector(aG,NSPDK);
		vector<unsigned> first_endpoint_list= vector<unsigned>();
		first_endpoint_list.push_back(i);
			//GetFirstEndpoints(aG, first_endpoint_list);

		int horizon=max(mDistance, mRadiusTwo); //   max(mDistance, mRadius);
			aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
			if (aG.Check() == false)
				throw logic_error("ERROR10: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

			InitFeatureCache(aG, mRadiusTwo );
			for (unsigned r = 0; r <= mRadius ; r++) { //max(1.0,mRadius / 2.0)
				for (unsigned d = 0; d <= mDistance; d++) { // era 0
					for (unsigned r2 = 0; r2 <= mRadiusTwo ; r2++) {
					SVector z;
					for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
						unsigned src_id = first_endpoint_list[i];
						if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
							SVector zv;

							vector<unsigned> endpoint_list(5);
							endpoint_list[0] = r2;
							endpoint_list[1] = d;

						//	unsigned src_code = GenerateVertexNeighbourhoodHashCode(aSrcID, aG, aRadius);
				//	TODO	for(unsigned src_code in DD)

							vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(i, d);
							for (unsigned dest_j = 0; dest_j < dest_id_list.size(); dest_j++) {
								unsigned dest_id = dest_id_list[dest_j];
								unsigned dest_code = 0;
								if (aG.GetVertexKernelPoint(dest_id) && aG.GetVertexAlive(dest_id)) { //proceed to extract features only if the *dest* vertex is a kernel point and is alive
									dest_code = GenerateVertexNeighbourhoodHashCode(dest_id, aG, r2); //
									// first is the DDK feature, second NSPDK feature
						//  for(const SVector::Pair *p = DD; p->i>=0; p++) {
									for(int distance=0; distance<features[r].size(); distance++){
										endpoint_list[4]=distance; //from root in DD
						 for(const SVector::Pair *p = features[r][distance]; p->i>=0; p++) {

									endpoint_list[2] = p->i;
										//TODO scorrere le feature di DD
										endpoint_list[3] = dest_code;
										//TEST output
								//		cout<<"R "<<endpoint_list[0]<<" d "<<endpoint_list[1]<<" DD "<<endpoint_list[2]<<" NSPDK "<<endpoint_list[3]<<" dRoot "<<endpoint_list[4]<<endl;
									unsigned code = HashFunc(endpoint_list, mHashBitMask);

									SVector z;
									z.set(code, p->v);
									zv.add(z);
								}
								}

								}
							}



							z.add(zv);
						}
						}

					if (mNormalization)
						z.normalize();
				out->add(z);
				} //d
			} //r
			} //r2







	//	cout<<"features size DD:"<<DD.sparse_size()<<" out:"<<out->sparse_size()<<endl;
		features.clear();


	return *out;
}

void ANSDDkernel_FeatureGenerator::GenerateAbstractVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x) {
	//cout<<"ABSTRACT FEATURES"<<endl;
	vector<unsigned> endpoint_list(5);
	for (unsigned r = 0; r <= mRadius ; r++) { //max(1.0,mRadius / 2.0)
					for (unsigned d = 0; d <= mDistance; d++) { // era 0
						for (unsigned r2 = 0; r2 <= mRadiusTwo ; r2++) { //max(1.0,mRadius / 2.0)


	endpoint_list[0] = r;
	endpoint_list[1] = d;

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
	//----------------------------

	//----------------------------------------
						//starting from each upper level vertex find all vertices at distance aDistance
	if (abstraction_of_list.size() == 0)
		throw std::logic_error("ERROR5: Something went wrong: expecting a non empty abstraction_of list.");
	set<unsigned> abstraction_of_set;
	for (unsigned i = 0; i < abstraction_of_list.size(); ++i) {
		unsigned v = abstraction_of_list[i];
		vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(v, d);
		//-----------?????????????????????????????
	//	abstraction_of_set.insert(abstraction_of_list.begin(), abstraction_of_list.end());
		abstraction_of_set.insert(dest_id_list.begin(), dest_id_list.end());

	}
	//starting from each lower level vertex find all vertices at distance aDistance
	if (part_of_list.size() == 0)
		throw std::logic_error("ERROR4: Something went wrong: expecting a non empty part_of list.");
	set<unsigned> part_of_set;
	for (unsigned i = 0; i < part_of_list.size(); ++i) {
		unsigned v = part_of_list[i];
		vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(v, d);
		//-----------?????????????????????????????
		part_of_set.insert(dest_id_list.begin(), dest_id_list.end());
	}
	//make all possible pairs of distant lower level vertices with distant upper level vertices
	for (set<unsigned>::iterator it = part_of_set.begin(); it != part_of_set.end(); ++it) {
		unsigned part_of_id = *it;

		//TODO modifica con mio kernel
		vector< vector< SVector> > features; //(dag_h+1);
		SVector DD= generate_vertex_feature_vector_DIVIDED_DISTANCE(part_of_id, aG ,aRadius,features);
		//SVector DD= generate_vertex_feature_vector_DIVIDED(part_of_id, aG ,aRadius,features);

		//For every feature in DD, by height and distance from the root
		for(int distance=0; distance<features[r].size(); distance++){
			endpoint_list[4]=distance; //from root in DD
			for(const SVector::Pair *p = features[r][distance]; p->i>=0; p++) {

			endpoint_list[2] = p->i;
			//TODO scorrere le feature di DD
				//TEST output
			//		cout<<"R "<<endpoint_list[0]<<" d "<<endpoint_list[1]<<" DD "<<endpoint_list[2]<<" NSPDK "<<endpoint_list[3]<<" dRoot "<<endpoint_list[4]<<endl;

			//
		for (set<unsigned>::iterator jt = abstraction_of_set.begin(); jt != abstraction_of_set.end(); ++jt) {
			unsigned abstraction_of_id = *jt;
			//build features with one neighborhood graph signature from the lower level and one from the upper
			unsigned abstraction_of_code = GenerateVertexNeighbourhoodHashCode(abstraction_of_id, aG, r2);
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
	}


	//i have to repeat for both sides for ddk
	//make all possible pairs of distant lower level vertices with distant upper level vertices
	for (set<unsigned>::iterator it = part_of_set.begin(); it != part_of_set.end(); ++it) {
		unsigned part_of_id = *it;

		//TODO modifica con mio kernel
		vector< vector< SVector> > features; //(dag_h+1);
		SVector DD= generate_vertex_feature_vector_DIVIDED_DISTANCE(part_of_id, aG ,aRadius,features);

			//For every feature in DD, by height and distance from the root
			for(int distance=0; distance<features[r].size(); distance++){
				endpoint_list[4]=distance; //from root in DD
				for(const SVector::Pair *p = features[r][distance]; p->i>=0; p++) {

				endpoint_list[2] = p->i;

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

}
					}
	}
}
}



