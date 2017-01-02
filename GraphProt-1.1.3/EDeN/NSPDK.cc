#include "Utility.h"
#include "BaseGraphClass.h"
#include "GraphClass.h"
#include "NSPDK_FeatureGenerator.h"
#include "wrapper.h"
#include "vectors.h"
#include <ctime>
#include <numeric>
#include <list>
#include <stdio.h>

using namespace std;

//---------------------------------------------------------------------------------
string DATE(__DATE__);
const string CREDITS("Name: NeighborhoodSubgraphPaiwiseDistanceKernel\nVersion: 9.2\nProgrammer: Fabrizio Costa\nDate:3 July 2012");
//---------------------------------------------------------------------------------
ofstream LOGF("log", std::ios_base::app);

FlagsService& The_FlagsService = FlagsService::get_instance();

class ParameterWrapperClass {
public:
	ParameterWrapperClass() :
			mGspanInputFileName(""), mSparseBinaryInputFileName(""), mSparseASCIIInputFileName(""), mSparseASCIITestInputFileName(""), mTrainTargetInputFileName(""), mRadiusMax(1), mDistanceMax(4), mMatchingType("hard"), mVerbose(false), mFeatureBitSize(30), mMinKernel(false), mType("nspdk"), mGraphType("UNDIRECTED"), mNormalization(true), mNumNearestNeighbors(10), mNumHashFunctions(250), mSampleSize(10), mNonRedundantFilter(1), mOutputFeatures(false), mOutputAccuracy(false), mOutputFeatureMap(false), mOutputKernel(false), mOutputApproximateKNN(false), mOutputTrueKNN(false), mOutputCluster(false), mOutputApproximateCluster(false), mOutputApproximateKNNPrediction(false), mOutputTrueKNNPrediction(false), mOutputHashEncoding(false), mCaching(true), mTrueSort(true), mEccessNeighbourSizeFactor(10), mHashFactor(1), mNumCenters(10), mSizeThreshold(30), mImbalanceTolerance(1.5), mWhiteListFileName(""), mBlackListFileName(""), mGreyListFileName(""), mMaxIntersectionSize(0), mFractionCenterScan(1), mNumMinHashFunc(50), mMaxSizeBin(.1), mRandSeed(1), mSharedNeighbourhood(false), mFullDensityEstimation(false), mMaxRefinement(100000), mDebug(0) {
	}

	void Usage(string aCommandName) {
		cerr
				<< "Usage: "
				<< aCommandName
				<< endl
				<< "-fg <file name gspan format> for input from file "
				<< endl
				<< "-fsb <file name sparse format binary> for input from file"
				<< endl
				<< "-fsa <file name sparse format ascii> for input from file"
				<< endl
				<< "-fsats <file name sparse format ascii> for input from file for test set"
				<< endl
				<< "-ftrt <file name> for target for train set"
				<< endl
				<< "[-knn <num nearest neighbors> (default: "
				<< mNumNearestNeighbors
				<< ")]"
				<< endl
				<< "-wl <file name for white listed instances ids (1 based)> for data set"
				<< endl
				<< "-bl <file name for black listed instances ids (1 based)> for data set"
				<< endl
				<< "-gl <file name for grey listed instances ids (1 based)> for data set"
				<< endl
				<< "[-oacc flag to output accuracy of approximation (default: "
				<< mOutputAccuracy
				<< ")]"
				<< endl
				<< "[-of flag to output feature encoding (default: "
				<< mOutputFeatures
				<< ")]"
				<< endl
				<< "[-ofm flag to output feature map encoding (default: "
				<< mOutputFeatureMap
				<< ")]"
				<< endl
				<< "[-ok flag to output kernel matrix (default: "
				<< mOutputKernel
				<< ")]"
				<< endl
				<< "[-oaknn flag to output approximate k-nearest neighburs (default: "
				<< mOutputApproximateKNN
				<< ")]"
				<< endl
				<< "[-otknn flag to output true (i.e. implies full kernel matrix evaluation) k-nearest neighburs (default: "
				<< mOutputTrueKNN
				<< ")]"
				<< endl
				<< "[-oc flag to output clusters (default: "
				<< mOutputCluster
				<< ")]"
				<< endl
				<< "[-oac flag to output approximate clusters (default: "
				<< mOutputApproximateCluster
				<< ")]"
				<< endl
				<< "[-nc num centers (default: "
				<< mNumCenters
				<< ")]"
				<< endl
				<< "[-fcs fraction of dataset size to scan for centers (default: "
				<< mFractionCenterScan
				<< ")]"
				<< endl
				<< "[-nhf <num hash functions> for the Locality Sensitive Hashing function (default: "
				<< mNumHashFunctions
				<< ")]"
				<< endl
				<< "[-hf <hash factor> number of signatures to collate (default: "
				<< mHashFactor
				<< ")]"
				<< endl
				<< "[-msb <max size bin > (default: "
				<< mMaxSizeBin
				<< ") (expressed as a fraction of the dataset size)] "
				<< endl
				<< "[-ensf <eccess neighbour size factor> (default: "
				<< mEccessNeighbourSizeFactor
				<< ") (0 to avoid trimming)] "
				<< endl
				<< "[-ss <sample size> for clustering procedure (default: "
				<< mSampleSize
				<< ")]"
				<< endl
				<< "[-mi <maximum number of elements in intersection between center neighborhoods> (default: "
				<< mMaxIntersectionSize
				<< ")]"
				<< endl
				<< "[-oaknnp flag to output approximate knn prediction (default: "
				<< mOutputApproximateKNNPrediction
				<< ")]"
				<< endl
				<< "[-otknnp flag to output true knn prediction (default: "
				<< mOutputTrueKNNPrediction
				<< ")]"
				<< endl
				<< "[-ohe flag to output hash encoding (default: "
				<< mOutputHashEncoding
				<< ")]"
				<< endl
				<< "[-b <feature space bits size> (default: "
				<< mFeatureBitSize
				<< ")]"
				<< endl
				<< "[-R <max radius> (default: "
				<< mRadiusMax
				<< ")]"
				<< endl
				<< "[-D <max distance relations> (default: "
				<< mDistanceMax
				<< ")]"
				<< endl
				<< "[-gt <graph type DIRECTED|UNDIRECTED> (default: "
				<< mGraphType
				<< ")]"
				<< endl
				<< "[-anhf <number of hash functions for abstract> (default: "
				<< mNumMinHashFunc
				<< ")]"
				<< endl
				<< "[-T <nspdk> (default: "
				<< mType
				<< ")]"
				<< endl
				<< "[-t <hard | soft | hard_soft> as neighborhood matching use HARD=exact matching, use SOFT=attribute matching with root identifier as radius 0 neighborhood, use HARD-SOFT=attribute matching with root identifier as full neighborhood encoding (default: "
				<< mMatchingType
				<< ")]"
				<< endl
				<< "[-mink flag to set minimum kernel rather than dot product (default: "
				<< mMinKernel
				<< ")]"
				<< endl
				<< "[-nn flag to de-acivate normalization (default: "
				<< !mNormalization
				<< ")]"
				<< endl
				<< "[-nrt <similarity filtering of redundant centers [0,1]> for clustering procedure (default: "
				<< mNonRedundantFilter
				<< ") (the smaller the less similar the centers)]"
				<< endl
				<< "[-no-cache flag to deactivate caching of kernel value computation (to minimize memory usage) (default: "
				<< !mCaching
				<< ")]"
				<< endl
				<< "[-no-true-sort flag to deactivate sorting approximate neighbours with true kernel computation (default: "
				<< !mTrueSort
				<< ")]"
				<< endl
				<< "[-st <size threshold> (default: "
				<< mSizeThreshold
				<< ")]"
				<< endl
				<< "[-it <imbalance tolerance> (default: "
				<< mImbalanceTolerance
				<< ")]"
				<< endl
				<< "[-v flag for verbose output (default: "
				<< mVerbose
				<< ")]"
				<< endl
				<< "[-rs rand seed (default: "
				<< mRandSeed
				<< ")]"
				<< endl
				<< "[-usn use shared neighbourhood to weight center density (default: "
				<< mSharedNeighbourhood
				<< ")]"
				<< endl
				<< "[-fde perform full density estimation (default: "
				<< mFullDensityEstimation
				<< ")]"
				<< endl
				<< "[-mri <num max refinement iterations> (default: "
				<< mMaxRefinement
				<< ")]"
				<< endl
				<< "[-debug <debug level> for NSPDK data structures (default: "
				<< mDebug
				<< ")]"
				<< endl;
		exit(0);
	}

	void Init(int argc, char** argv) {
		vector<string> options;
		for (int i = 1; i < argc; i++)
			options.push_back(argv[i]);
		for (vector<string>::iterator it = options.begin(); it != options.end(); ++it) {
			if ((*it) == "-h" || (*it) == "--help") Usage(argv[0]);
			else if ((*it) == "-fg") mGspanInputFileName = (*(++it));
			else if ((*it) == "-fsb") mSparseBinaryInputFileName = (*(++it));
			else if ((*it) == "-fsa") mSparseASCIIInputFileName = (*(++it));
			else if ((*it) == "-fsats") mSparseASCIITestInputFileName = (*(++it));
			else if ((*it) == "-ftrt") mTrainTargetInputFileName = (*(++it));
			else if ((*it) == "-wl") mWhiteListFileName = (*(++it));
			else if ((*it) == "-bl") mBlackListFileName = (*(++it));
			else if ((*it) == "-gl") mGreyListFileName = (*(++it));
			else if ((*it) == "-knn") mNumNearestNeighbors = stream_cast<unsigned>(*(++it));
			else if ((*it) == "-ensf") mEccessNeighbourSizeFactor = stream_cast<double>(*(++it));
			else if ((*it) == "-hf") mHashFactor = stream_cast<unsigned>(*(++it));
			else if ((*it) == "-R") mRadiusMax = stream_cast<double>(*(++it));
			else if ((*it) == "-D") mDistanceMax = stream_cast<double>(*(++it));
			else if ((*it) == "-t") mMatchingType = (*(++it));
			else if ((*it) == "-v") mVerbose = true;
			else if ((*it) == "-mink") mMinKernel = true;
			else if ((*it) == "-b") mFeatureBitSize = stream_cast<int>(*(++it));
			else if ((*it) == "-T") mType = (*(++it));
			else if ((*it) == "-gt") mGraphType = (*(++it));
			else if ((*it) == "-of") mOutputFeatures = true;
			else if ((*it) == "-oacc") mOutputAccuracy = true;
			else if ((*it) == "-ofm") mOutputFeatureMap = true;
			else if ((*it) == "-ok") mOutputKernel = true;
			else if ((*it) == "-oaknn") mOutputApproximateKNN = true;
			else if ((*it) == "-otknn") mOutputTrueKNN = true;
			else if ((*it) == "-oaknnp") mOutputApproximateKNNPrediction = true;
			else if ((*it) == "-otknnp") mOutputTrueKNNPrediction = true;
			else if ((*it) == "-ohe") mOutputHashEncoding = true;
			else if ((*it) == "-oc") mOutputCluster = true;
			else if ((*it) == "-oac") mOutputApproximateCluster = true;
			else if ((*it) == "-nn") mNormalization = false;
			else if ((*it) == "-debug") mDebug = stream_cast<int>(*(++it));
			else if ((*it) == "-nhf") mNumHashFunctions = stream_cast<unsigned>(*(++it));
			else if ((*it) == "-ss") mSampleSize = stream_cast<unsigned>(*(++it));
			else if ((*it) == "-nrt") mNonRedundantFilter = stream_cast<double>(*(++it));
			else if ((*it) == "-no-cache") mCaching = false;
			else if ((*it) == "-no-true-sort") mTrueSort = false;
			else if ((*it) == "-nc") mNumCenters = stream_cast<unsigned>(*(++it));
			else if ((*it) == "-fcs") mFractionCenterScan = stream_cast<double>(*(++it));
			else if ((*it) == "-st") mSizeThreshold = stream_cast<unsigned>(*(++it));
			else if ((*it) == "-it") mImbalanceTolerance = stream_cast<double>(*(++it));
			else if ((*it) == "-mi") mMaxIntersectionSize = stream_cast<unsigned>(*(++it));
			else if ((*it) == "-anhf") mNumMinHashFunc = stream_cast<unsigned>(*(++it));
			else if ((*it) == "-msb") mMaxSizeBin = stream_cast<double>(*(++it));
			else if ((*it) == "-rs") mRandSeed = stream_cast<unsigned>(*(++it));
			else if ((*it) == "-usn") mSharedNeighbourhood = true;
			else if ((*it) == "-fde") mFullDensityEstimation = true;
			else if ((*it) == "-mr") mMaxRefinement = stream_cast<unsigned>(*(++it));
			else {
				cerr << "Unrecognized parameter: " << (*it) << "." << endl;
				throw exception();
			}
		}
		if (!(mMatchingType == "hard" || mMatchingType == "soft" || mMatchingType == "hard_soft" || mMatchingType == "multiview" || mMatchingType == "mixed")) {
			cerr << "Wrong value for parameter: -t: " << mMatchingType << endl;
			throw exception();
		}
	}

public:
	string mGspanInputFileName;
	string mSparseBinaryInputFileName;
	string mSparseASCIIInputFileName;
	string mSparseASCIITestInputFileName;
	string mTrainTargetInputFileName;
	double mRadiusMax;
	double mDistanceMax;
	string mMatchingType;
	bool mVerbose;
	int mFeatureBitSize;
	bool mMinKernel;
	string mType;
	string mGraphType;
	bool mNormalization;
	unsigned mNumNearestNeighbors;
	unsigned mNumHashFunctions;
	unsigned mSampleSize;
	double mNonRedundantFilter;
	bool mOutputFeatures;
	bool mOutputAccuracy;
	bool mOutputFeatureMap;
	bool mOutputKernel;
	bool mOutputApproximateKNN;
	bool mOutputTrueKNN;
	bool mOutputCluster;
	bool mOutputApproximateCluster;
	bool mOutputApproximateKNNPrediction;
	bool mOutputTrueKNNPrediction;
	bool mOutputHashEncoding;
	bool mCaching;
	bool mTrueSort;
	double mEccessNeighbourSizeFactor;
	unsigned mHashFactor;
	unsigned mNumCenters;
	unsigned mSizeThreshold;
	double mImbalanceTolerance;
	string mWhiteListFileName;
	string mBlackListFileName;
	string mGreyListFileName;
	unsigned mMaxIntersectionSize;
	double mFractionCenterScan;
	unsigned mNumMinHashFunc;
	double mMaxSizeBin;
	unsigned mRandSeed;
	bool mSharedNeighbourhood;
	bool mFullDensityEstimation;
	unsigned mMaxRefinement;
	int mDebug;
} PARAM_OBJ;

typedef std::tr1::unordered_map<unsigned, int> umap_uint_int;
typedef std::tr1::unordered_map<unsigned, vector<unsigned> > umap_uint_vec_uint;

class NSPDKClass {
protected:
	NSPDK_FeatureGenerator* pmFeatureGenerator;

	vector<SVector> mDataset;
	vector<umap_uint_vec_uint> mBinDataStructure;
	map<pair<unsigned, unsigned>, double> mKernelMap;
	//multimap<unsigned, unsigned> mInvertedIndex;
	umap_uint_vec_uint mSignatureMap;

	double mAlpha;
	vector<double> mApproximateDensityMap;
	vector<double> mTrueDensityMap;
	umap_uint_vec_uint mApproximateNeighborhoodMap;
	vector<unsigned> mIdMap;
	vector<bool> mFilteredHashFunctionList;
	vector<int> mGreyList;
public:
	NSPDKClass(NSPDK_FeatureGenerator* paFeatureGenerator) :
			pmFeatureGenerator(paFeatureGenerator) {
	}

	void Generate(const GraphClass& aG, SVector& oX) {
		//create base graph features
		if (PARAM_OBJ.mType == "abstnspdk") aG.ComputePairwiseDistanceInformation(PARAM_OBJ.mDistanceMax, PARAM_OBJ.mRadiusMax);
		pmFeatureGenerator->generate_feature_vector(aG, oX);
	}

	void InputStringList(const string& aFileName, vector<string>& oStringList) {
		cout << "Reading " << aFileName << endl;
		ifstream fin;
		fin.open(aFileName.c_str());
		if (!fin) throw range_error("Cannot open file:" + aFileName);
		ProgressBar progress_bar;
		while (!fin.eof() && fin.good()) {
			string line;
			getline(fin, line);
			stringstream ss;
			ss << line << endl;
			while (!ss.eof() && ss.good()) {
				string target;
				ss >> target;
				if (target != "") {
					oStringList.push_back(target);
					progress_bar.Count();
				}
			}
		}
		fin.close();
	}

	void InputIntList(const string& aFileName, vector<int>& oList) {
		cout << "Reading " << aFileName << endl;
		ifstream fin;
		fin.open(aFileName.c_str());
		if (!fin) throw range_error("Cannot open file:" + aFileName);
		ProgressBar progress_bar;
		while (!fin.eof() && fin.good()) {
			string line;
			getline(fin, line);
			if (line != "") {
				stringstream ss;
				ss << line;
				int target;
				ss >> target;
				oList.push_back(target);
				progress_bar.Count();
			}
		}
		fin.close();
	}

	void Input() {
		if (PARAM_OBJ.mSparseBinaryInputFileName != "") InputSparse(PARAM_OBJ.mSparseBinaryInputFileName, "binary");
		else if (PARAM_OBJ.mSparseASCIIInputFileName != "") InputSparse(PARAM_OBJ.mSparseASCIIInputFileName, "ascii");
		else if (PARAM_OBJ.mGspanInputFileName != "") {
			if (PARAM_OBJ.mOutputFeatures) {
				Load(PARAM_OBJ.mGspanInputFileName, "direct");
			} else if (PARAM_OBJ.mOutputApproximateCluster) {
				Load(PARAM_OBJ.mGspanInputFileName, "approximate");
			} else {
				Load(PARAM_OBJ.mGspanInputFileName, "memory");
			}
		} else throw range_error("ERROR:No input file name specified");
	}

//DirectProcess true=1 false=0 discard vector=2
	void Load(const string& aInputFileName, const string& aTypeOfProcess) {
		ofstream ofs_f;
		ofstream ofs_fb;
		if (aTypeOfProcess == "direct") {
			string ofname = aInputFileName + ".feature";
			ofs_f.open(ofname.c_str());
			ofname = aInputFileName + ".feature_bin";
			ofs_fb.open(ofname.c_str());
		}

		//read white list
		vector<int> select_list;
		if (PARAM_OBJ.mWhiteListFileName != "") {
			InputIntList(PARAM_OBJ.mWhiteListFileName, select_list);
		}
		//read black list
		if (PARAM_OBJ.mBlackListFileName != "") {
			InputIntList(PARAM_OBJ.mBlackListFileName, select_list);
		}
		//read grey list
		if (PARAM_OBJ.mGreyListFileName != "") {
			InputIntList(PARAM_OBJ.mGreyListFileName, mGreyList);
		}
		set<int> select_list_set;
		select_list_set.insert(select_list.begin(), select_list.end());

		cout << "Reading gspan data and computing features" << endl;
		ifstream fin;
		fin.open(aInputFileName.c_str());
		if (!fin) throw range_error("Cannot open file:" + aInputFileName);
		ProgressBar progress_bar;
		int counter = 1;
		while (!fin.eof()) {
			GraphClass G;
			SetGraphFromFileGSPAN(fin, G);
			SVector x;

			//only if counter id is consistent with white and black list (if they have been specified) then accept the instance
			bool accept_flag = true;
			if (PARAM_OBJ.mWhiteListFileName != "") {
				if (select_list_set.count(counter) > 0) accept_flag = true;
				else accept_flag = false;
			}
			if (PARAM_OBJ.mBlackListFileName != "") {
				if (select_list_set.count(counter) > 0) accept_flag = false;
				else accept_flag = true;
			}

			if (accept_flag == true) {
				Generate(G, x);
				if (aTypeOfProcess == "direct") {
					ofs_f << x;
					x.save(ofs_fb);
				} else if (aTypeOfProcess == "approximate") {
					AddToBinDataStructure(x);
				} else if (aTypeOfProcess == "memory") {
					mDataset.push_back(x);
				} else throw range_error("ERROR:Invalid load mode: <" + aTypeOfProcess + ">");
				mIdMap.push_back(counter);
			}
			progress_bar.Count();
			counter++;
		}
		fin.close();
	}

	void InputSparse(const string& aInputFileName, string aMode) {
		ifstream fin;
		fin.open(aInputFileName.c_str());
		if (!fin) throw range_error("Cannot open file:" + aInputFileName);
		InputSparse(fin, aMode, mDataset);
		fin.close();
	}

	void InputSparse(const string& aInputFileName, string aMode, vector<SVector>& oDataset) {
		ifstream fin;
		fin.open(aInputFileName.c_str());
		if (!fin) throw range_error("Cannot open file:" + aInputFileName);
		InputSparse(fin, aMode, oDataset);
		fin.close();
	}

	void InputSparse(ifstream& aFin, string aMode, vector<SVector>& oDataset) {
		//read white list
		vector<int> select_list;
		if (PARAM_OBJ.mWhiteListFileName != "") {
			InputIntList(PARAM_OBJ.mWhiteListFileName, select_list);
		}
		//read black list
		if (PARAM_OBJ.mBlackListFileName != "") {
			InputIntList(PARAM_OBJ.mBlackListFileName, select_list);
		}
		//read grey list
		if (PARAM_OBJ.mGreyListFileName != "") {
			InputIntList(PARAM_OBJ.mGreyListFileName, mGreyList);
		}
		set<int> select_list_set;
		select_list_set.insert(select_list.begin(), select_list.end());

		cout << "Reading file in " << aMode << " mode" << endl;
		int counter = 1;
		ProgressBar progress_bar;
		while (!aFin.eof() && aFin.good()) {
			SVector x;
			if (aMode == "binary") x.load(aFin);
			else ParseASCIILine2Vector(aFin, x);
			if (InstanceIsValid(x) == true) {
				//only if counter id is consistent with white and black list (if they have been specified) then accept the instance
				bool accept_flag = true;
				if (PARAM_OBJ.mWhiteListFileName != "") {
					if (select_list_set.count(counter) > 0) accept_flag = true;
					else accept_flag = false;
				}
				if (PARAM_OBJ.mBlackListFileName != "") {
					if (select_list_set.count(counter) > 0) accept_flag = false;
					else accept_flag = true;
				}
				if (accept_flag == true) {
					if (PARAM_OBJ.mOutputApproximateCluster) {
						AddToBinDataStructure(x);
					} else {
						oDataset.push_back(x);
					}
					mIdMap.push_back(counter);
				}
				progress_bar.Count();
				counter++;
			} else {
			} //discard non valid instances
		}
	}

	inline void ParseASCIILine2Vector(ifstream& aFin, SVector& aX) {
		string line;
		getline(aFin, line);
		if (line == "") return;
		stringstream ss;
		ss << line << endl;
		while (!ss.eof() && ss.good()) {
			string key_value;
			ss >> key_value;
			size_t limit = key_value.find_first_of(":", 0);
			if (limit != string::npos) { //if the delimiter ':' is found then proceed
				string key = key_value.substr(0, limit);
				string value = key_value.substr(limit + 1, key_value.size());
				unsigned key_int = stream_cast<unsigned>(key);
				double val_real = stream_cast<double>(value);
				aX.set(key_int, val_real);
			}
		}
	}

	inline bool InstanceIsValid(SVector& aX) {
		bool is_valid = false;
		//if (dot(aX,aX)>0) is_valid=true;
		if (aX.sparse_size() > 0) is_valid = true;
		return is_valid;
	}

	void SetGraphFromFileGSPAN(istream& in, GraphClass& oG) {
		//status
		vector<bool> vertex_status;
		vertex_status.push_back(true); //kernel point
		vertex_status.push_back(true); //kind
		vertex_status.push_back(true); //viewpoint
		vertex_status.push_back(false); //dead
		vertex_status.push_back(false); //abstraction

		vector<bool> edge_status;
		edge_status.push_back(false); //edge dead
		edge_status.push_back(false); //edge abstraction_of
		edge_status.push_back(false); //edge part_of

		map<string, int> index_map_nominal_to_real;
		string line;
		getline(in, line);
		assert(line[0]=='t');
		//first line must have as first char a 't'
		static unsigned line_count = 1; //line counter with overall file scope
		while (!in.eof() && in.good() && in.peek() != 't' && getline(in, line)) { //read until next 't' or end of file
			line_count++;
			stringstream ss;
			ss << line << endl;
			char code;
			ss >> code;
			if (code == 'v') {
				//extract vertex id and make map nominal_id -> real_id
				string nominal_vertex_index;
				ss >> nominal_vertex_index;
				unsigned real_vertex_index = oG.InsertVertex();
				index_map_nominal_to_real[nominal_vertex_index] = real_vertex_index;
				//label
				vector<string> vertex_symbolic_attribute_list;
				string label;
				ss >> label;
				vertex_symbolic_attribute_list.push_back(label);

				oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
				oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);
				char vertex_abstraction_code = label.at(0);
				if (vertex_abstraction_code == '^') {
					oG.SetVertexKernelPoint(real_vertex_index, false);
					oG.SetVertexAbstraction(real_vertex_index, true);
				}

			} else if (code == 'e') {
				//extract src and dest vertex id
				string nominal_src_index, nominal_dest_index;
				string label;
				ss >> nominal_src_index >> nominal_dest_index >> label;
				if (index_map_nominal_to_real.count(nominal_src_index) == 0) throw range_error("Error in line:" + stream_cast<string>(line_count) + " What: Edge with source endpoint in non decleared vertex with id " + nominal_src_index);
				if (index_map_nominal_to_real.count(nominal_dest_index) == 0) throw range_error("Error in line:" + stream_cast<string>(line_count) + " What: Edge with destination endpoint in non decleared vertex with id " + nominal_dest_index);
				vector<string> edge_symbolic_attribute_list;
				edge_symbolic_attribute_list.push_back(label);
				unsigned real_src_index = index_map_nominal_to_real[nominal_src_index];
				unsigned real_dest_index = index_map_nominal_to_real[nominal_dest_index];
				unsigned edge_index = oG.InsertEdge(real_src_index, real_dest_index);
				oG.SetEdgeSymbolicAttributeList(edge_index, edge_symbolic_attribute_list);

				vector<bool> current_edge_status(edge_status);
				char edge_abstraction_code = label.at(0);
				if (edge_abstraction_code == '^') current_edge_status[1] = true;
				if (edge_abstraction_code == '@') current_edge_status[2] = true;
				oG.SetEdgeStatusAttributeList(edge_index, current_edge_status);

				if (PARAM_OBJ.mGraphType == "UNDIRECTED") {
					unsigned reverse_edge_index = oG.InsertEdge(real_dest_index, real_src_index);
					oG.SetEdgeSymbolicAttributeList(reverse_edge_index, edge_symbolic_attribute_list);
					oG.SetEdgeStatusAttributeList(reverse_edge_index, current_edge_status);
				}
			} else {
			} //NOTE: ignore other markers
		}
		if (PARAM_OBJ.mGraphType == "DIRECTED") {
			unsigned vsize = oG.VertexSize();
			//add a copy of all vertices
			for (unsigned i = 0; i < vsize; i++) {
				unsigned real_vertex_index = oG.InsertVertex();
				assert(real_vertex_index=i+vsize);
				vector<string> r_vertex_symbolic_attribute_list = oG.GetVertexSymbolicAttributeList(i);
				for (unsigned t = 0; t < r_vertex_symbolic_attribute_list.size(); t++) //prepend a prefix to mark the reverse direction
					r_vertex_symbolic_attribute_list[t] = "r." + r_vertex_symbolic_attribute_list[t];
				oG.SetVertexSymbolicAttributeList(real_vertex_index, r_vertex_symbolic_attribute_list);
				oG.SetVertexStatusAttributeList(real_vertex_index, oG.GetVertexStatusAttributeList(i)); //assign original status vector
			}
			//copy all edges swapping src with dest
			for (unsigned i = 0; i < vsize; i++) {
				//get all edges
				vector<unsigned> adj = oG.GetVertexAdjacentList(i);
				for (unsigned j = 0; j < adj.size(); j++) {
					unsigned orig_src = i;
					unsigned orig_dest = adj[j];
					unsigned reverse_src = orig_dest + vsize;
					unsigned reverse_dest = orig_src + vsize;
					unsigned edge_index = oG.InsertEdge(reverse_src, reverse_dest);
					oG.SetEdgeSymbolicAttributeList(edge_index, oG.GetEdgeSymbolicAttributeList(orig_src, orig_dest));
					oG.SetEdgeStatusAttributeList(edge_index, oG.GetEdgeStatusAttributeList(orig_src, orig_dest));
				}
			}
		}

	}

	void ComputeBinDataStructure() {
		string ofname = "hash_encoding";
		ofstream of(ofname.c_str());

		cout << "Computing bin data structure..." << endl;
		ProgressBar progress_bar;

		InitBinDataStructure();

		//fill structure
		for (unsigned i = 0; i < mDataset.size(); ++i) {
			vector<unsigned> min_list = ComputeHashSignature(i);
			if (PARAM_OBJ.mOutputHashEncoding) {
				for (unsigned j = 0; j < min_list.size(); j++)
					of << min_list[j] << " ";
				of << endl;
			}
			for (unsigned k = 0; k < PARAM_OBJ.mNumHashFunctions; ++k) {
				if (mBinDataStructure[k].count(min_list[k]) > 0) {
					mBinDataStructure[k][min_list[k]].push_back(i);
				} else {
					vector<unsigned> tmp;
					tmp.push_back(i);
					mBinDataStructure[k].insert(make_pair(min_list[k], tmp));
				}
			}
			//	mBinDataStructure[k].insert(make_pair(min_list[k], i));
			progress_bar.Count();
		}
	}

	void InitBinDataStructure() {
		//init structure
		mBinDataStructure.clear();
		for (unsigned k = 0; k < PARAM_OBJ.mNumHashFunctions; ++k)
			mBinDataStructure.push_back(umap_uint_vec_uint());
	}

	void AddToBinDataStructure(SVector& aX) {
		static unsigned id_counter = 0;
		vector<unsigned> min_list = ComputeHashSignature(aX, id_counter);
		for (unsigned k = 0; k < PARAM_OBJ.mNumHashFunctions; ++k) {
			if (mBinDataStructure[k].count(min_list[k]) > 0) {
				mBinDataStructure[k][min_list[k]].push_back(id_counter);
			} else {
				vector<unsigned> tmp;
				tmp.push_back(id_counter);
				mBinDataStructure[k].insert(make_pair(min_list[k], tmp));
			}
		}
		id_counter++;
		mDataset.push_back(SVector());
	}

	inline vector<unsigned> ComputeHashSignature(unsigned aID) {
		if (mSignatureMap.count(aID) > 0) return mSignatureMap[aID];
		else {
			vector<unsigned> signature = ComputeHashSignature(mDataset[aID], aID);
			mSignatureMap[aID] = signature;
			return signature;
		}
	}

	inline vector<unsigned> ComputeHashSignature(SVector& aX, unsigned aID) {
		unsigned effective_num_hash_functions = PARAM_OBJ.mNumHashFunctions * PARAM_OBJ.mHashFactor;
		const unsigned MAXUNSIGNED = 2 << 30;
		vector<unsigned> signature;
		//prepare a vector containing the k min values
		for (unsigned k = 0; k < effective_num_hash_functions; ++k)
			signature.push_back(MAXUNSIGNED);
		unsigned size = (unsigned) aX.sparse_size();
		//for each element of the sparse vector
		for (unsigned f = 0; f < size; ++f) {
			//extract only the feature id (i.e. ignore the actial value)
			unsigned hash_id = aX.extract_component(f).first;
			if (hash_id == 0) {
				//feature is should not be 0 as the subsequent rehashing can encounter problems
				cout << "In sequence with id: " << mIdMap[aID] << endl; /////
				cout << "Warning: Feature ID = 0. Feature ID  should be strictly > 0" << endl;
				hash_id = 1; //force collision between feature 0 and 1
			}
			for (unsigned k = 0; k < effective_num_hash_functions; ++k) { //for all k hashes
				unsigned new_hash = IntHash(hash_id, MAXUNSIGNED, k); //rehash the feature id with a procedure that is aware of the index k
				if (signature[k] > new_hash) signature[k] = new_hash; //keep the minimum value only
			}
		}
		//compact signature
		vector<unsigned> compact_signature;
		for (unsigned i = 0; i < PARAM_OBJ.mNumHashFunctions; ++i)
			compact_signature.push_back(0);
		for (unsigned i = 0; i < signature.size(); ++i) {
			//add up several signatures to compose the new signature
			compact_signature[i % PARAM_OBJ.mNumHashFunctions] += signature[i];
		}
		return compact_signature;
	}

	inline vector<unsigned> ComputeHashSignatureNatural(SVector& aX) {
		unsigned effective_num_hash_functions = PARAM_OBJ.mNumHashFunctions * PARAM_OBJ.mHashFactor;
		vector<unsigned> signature;
		unsigned size = (unsigned) aX.sparse_size();
		unsigned effective_size = size < effective_num_hash_functions ? size : effective_num_hash_functions;
		for (unsigned f = 0; f < effective_size; ++f) {
			//do not rehash; instead use simply the f-th feature as
			//signature assuming they are indeed randomly distributed
			signature.push_back(aX.extract_component(f).first);
		}

		//compact signature
		vector<unsigned> compact_signature;
		for (unsigned i = 0; i < signature.size(); i = i + PARAM_OBJ.mHashFactor) {
			unsigned new_hash = 0;
			for (unsigned j = 0; j < PARAM_OBJ.mHashFactor; j++)
				new_hash += signature[i + j];
			compact_signature.push_back(new_hash);
		}
		return compact_signature;
	}

	void OutputBinDataStructureStatistics() const {
		VectorClass bin_size_stats, bin_type_stats;

		for (unsigned k = 0; k < mBinDataStructure.size(); ++k) {
			unsigned bin_type_counter = 0;
			for (umap_uint_vec_uint::const_iterator it = mBinDataStructure[k].begin(); it != mBinDataStructure[k].end(); ++it) {
				//unsigned key = it->first;
				vector<unsigned> dat = it->second;
				unsigned bin_size = dat.size();
				//unsigned  bin_size = mBinDataStructure[k][key].size();
				bin_size_stats.PushBack(bin_size);
				bin_type_counter++;
			}
			bin_type_stats.PushBack(bin_type_counter);
		}

		cout << "Num bins statistics: ";
		bin_type_stats.OutputStatistics(cout);
		cout << endl;
		cout << "Size bins statistics: ";
		bin_size_stats.OutputStatistics(cout);
		cout << endl;
	}

	vector<unsigned> ComputeApproximateNeighborhood(unsigned aID) {
		if (mApproximateNeighborhoodMap.count(aID) == 0) {
			vector<unsigned> hash_signature = ComputeHashSignature(aID);
			vector<unsigned> neighborhood = ComputeApproximateNeighborhood(hash_signature);
			//select neighborhood under true similarity function on the subset of indiced returned by ComputeApproximateNeighborhood
			vector<unsigned> true_neighborhood = ComputeTrueSubNeighborhood(aID, neighborhood);
			mApproximateNeighborhoodMap[aID] = true_neighborhood;
		}
		return mApproximateNeighborhoodMap[aID];
	}

	vector<unsigned> ComputeApproximateNeighborhood(const vector<unsigned>& aInstanceSignature) {
		umap_uint_int neighborhood;
		vector<pair<unsigned, double> > vec;
		for (unsigned k = 0; k < PARAM_OBJ.mNumHashFunctions; ++k) {
			unsigned hash_id = aInstanceSignature[k];
			unsigned collision_size = mBinDataStructure[k][hash_id].size();

			if (collision_size < PARAM_OBJ.mMaxSizeBin * mDataset.size()) {
				//fill neighborhood set counting number of occurrences
				for (vector<unsigned>::iterator it = mBinDataStructure[k][hash_id].begin(); it != mBinDataStructure[k][hash_id].end(); ++it) {
					unsigned instance_id = *it;
					if (neighborhood.count(instance_id) > 0) neighborhood[instance_id]++;
					else neighborhood[instance_id] = 1;
				}
			} else {
				mBinDataStructure[k].erase(hash_id);
			}
		}
		return TrimNeighborhood(neighborhood);
	}

	vector<unsigned> ComputeTrueSubNeighborhood(unsigned aID, vector<unsigned>& aApproximateNeighborhoodList) {
		vector<unsigned> neighbor;
		if (PARAM_OBJ.mTrueSort) {
			vector<pair<double, unsigned> > rank_list;
			for (unsigned i = 0; i < aApproximateNeighborhoodList.size(); ++i) {
				unsigned id_neighbour = aApproximateNeighborhoodList[i];
				double k = Kernel(aID, id_neighbour);
				rank_list.push_back(make_pair(-k, id_neighbour));
			}
			unsigned effective_size = min((unsigned) rank_list.size(), PARAM_OBJ.mNumNearestNeighbors);
			sort(rank_list.begin(), rank_list.end());

			for (unsigned j = 0; j < effective_size; j++) {
				neighbor.push_back(rank_list[j].second);
			}
		} else {
			unsigned effective_size = min((unsigned) aApproximateNeighborhoodList.size(), PARAM_OBJ.mNumNearestNeighbors);
			for (unsigned j = 0; j < effective_size; j++)
				neighbor.push_back(aApproximateNeighborhoodList[j]);
		}
		return neighbor;
	}

	vector<unsigned> ComputeTrueNeighborhood(unsigned aID) {
		const SVector& x = mDataset[aID];
		return ComputeTrueNeighborhood(x);
	}

	vector<unsigned> ComputeTrueNeighborhood(const SVector& aX) {
		vector<pair<double, unsigned> > rank_list;
		for (unsigned i = 0; i < mDataset.size(); ++i) {
			double k = dot(aX, mDataset[i]);
			rank_list.push_back(make_pair(-k, i));
		}
		unsigned effective_size = min((unsigned) rank_list.size(), PARAM_OBJ.mNumNearestNeighbors);
		partial_sort(rank_list.begin(), rank_list.begin() + effective_size, rank_list.end());
		vector<unsigned> neighbor;
		for (unsigned j = 0; j < effective_size; j++)
			neighbor.push_back(rank_list[j].second);
		return neighbor;
	}

	vector<unsigned> TrimNeighborhood(umap_uint_int& aNeighborhood) {
		const int MIN_BINS_IN_COMMON = 2; //Minimum number of bins that two instances have to have in common in order to be considered similar
		//given a list of neighbours with an associated occurences count, return only a fraction of the highest count ones
		vector<unsigned> neighborhood_list;
		if (PARAM_OBJ.mEccessNeighbourSizeFactor > 0) {
			//sort by num occurences
			vector<pair<int, unsigned> > count_list;
			for (umap_uint_int::const_iterator it = aNeighborhood.begin(); it != aNeighborhood.end(); ++it) {
				unsigned id = it->first;
				int count = it->second;
				if (count >= MIN_BINS_IN_COMMON) //NOTE: consider instances that have at least MIN_BINS_IN_COMMON
					count_list.push_back(make_pair(-count, id)); //NOTE:-count to sort from highest to lowest
			}
			sort(count_list.begin(), count_list.end());
			unsigned effective_size = min((unsigned) count_list.size(), (unsigned) (PARAM_OBJ.mEccessNeighbourSizeFactor * PARAM_OBJ.mNumNearestNeighbors));
			for (unsigned i = 0; i < effective_size; ++i)
				neighborhood_list.push_back(count_list[i].second);
		} else { //if mEccessNeighbourSizeFactor==0 then just consider all the ids in the approximate neighborhood
			for (umap_uint_int::const_iterator it = aNeighborhood.begin(); it != aNeighborhood.end(); ++it) {
				neighborhood_list.push_back(it->first);
			}
		}
		return neighborhood_list;
	}

	double ComputeApproximateDensity(unsigned aID) {
		double density = 0;
		unsigned i = aID;
		vector<unsigned> approximate_neighborhood = ComputeApproximateNeighborhood(i);
		if (mApproximateDensityMap[i] == -1) {
			if (PARAM_OBJ.mFullDensityEstimation == true) {
				for (unsigned it = 0; it < approximate_neighborhood.size(); it++) {
					i = approximate_neighborhood[it];
					density += CoreComputeDensity(i, approximate_neighborhood);
				}
				density = density / (approximate_neighborhood.size());
			} else {
				density = CoreComputeDensity(i, approximate_neighborhood);
			}
			mApproximateDensityMap[i] = density;
		} else {
			density = mApproximateDensityMap[i];
		}
		return density;
	}

	double CoreComputeDensity(unsigned aID, vector<unsigned>& aApproximateNeighborhood) {
		double density = 0;
		//compute kernel pairs between i and all elements in aApproximateNeighborhood
		for (unsigned j = 0; j < aApproximateNeighborhood.size(); j++) {
			unsigned u = aID;
			unsigned v = aApproximateNeighborhood[j];
			if (u != v) {
				double k_uv = Kernel(u, v);
				if (PARAM_OBJ.mSharedNeighbourhood) //if we use the shared neighborhood weighting than we multiply the similarity by a corrective factor given by the fraction of shared neighbours between the two instances
					k_uv *= ComputeSharedNeighborhoodSimilarity(u, v);
				density += k_uv;
			}
		}
		density = density / (aApproximateNeighborhood.size() - 1);
		return density;
	}

	double ComputeTrueDensity(unsigned aID) {
		double density = 0;
		if (mTrueDensityMap[aID] == -1) {
			vector<pair<double, unsigned> > sim_list;
			for (unsigned j = 0; j < mDataset.size(); j++) {
				if (aID != j) {
					double k_ij = Kernel(aID, j);
					density += k_ij;
				}
			}
			density = density / (mDataset.size() - 1);
			mTrueDensityMap[aID] = density;
		} else density = mTrueDensityMap[aID];
		return density;
	}

	double ComputeAverageDensity(set<unsigned>& aSet) {
		double density = 0;
		for (set<unsigned>::iterator it = aSet.begin(); it != aSet.end(); ++it) {
			unsigned id = *it;
			density += ComputeApproximateDensity(id);
		}
		return density;
	}
	vector<unsigned> ComputeMinimallyOverlappingHighDensityCenterList(unsigned aSampleSize, double aFractionCenterScan, unsigned aMaxIntersectionSize) {
		//compute density estimate for random fraction of instances in dataset
		unsigned data_size = mGreyList.size() > 0 ? mGreyList.size() : mDataset.size();
		unsigned effective_size = floor(data_size * aFractionCenterScan);
		cout << "Computing approximate density information for random sample of " << effective_size << " instances [" << aFractionCenterScan * 100 << "% of total size " << data_size << "]." << endl; ////

		//random selection of instances
		vector<unsigned> index_list;

		//select the candidate centers from gray list or from all ids
		if (mGreyList.size() > 0) { //if gray list is available then use ids from there
			for (unsigned i = 0; i < data_size; ++i)
				index_list.push_back(mGreyList[i] - 1); //NOTE -1 to adjust for the base 1
		} else { //else use all ids in dataset
			for (unsigned i = 0; i < data_size; ++i)
				index_list.push_back(i);
		}

		vector<unsigned> selected_index_list;
		//select either all the available centers
		if (aFractionCenterScan == 1) {
			selected_index_list.insert(selected_index_list.begin(), index_list.begin(), index_list.end());
		} else {
			//or select a random subset of candidate centers of size effective_size
			for (unsigned i = 0; i < data_size; ++i) {
				unsigned j = randomUnsigned(data_size);
				swap(index_list[i], index_list[j]);
			}
			for (unsigned i = 0; i < effective_size; ++i)
				selected_index_list.push_back(index_list[i]);
		}

		//compute density estimate
		vector<pair<double, unsigned> > density_list;
		{
			ProgressBar progress_bar;
			for (unsigned j = 0; j < selected_index_list.size(); ++j) {
				unsigned i = selected_index_list[j];
				double density = ComputeApproximateDensity(i);
				density_list.push_back(make_pair(-density, i));
				progress_bar.Count();

			}
		}

		//select non overlapping centers in decreasing order of density
		sort(density_list.begin(), density_list.end());
		vector<unsigned> result;
		set<unsigned> active_neighborhood;
		{
			ProgressBar progress_bar;
			cout << "Computing minimally overlapping high density center list for up to " << aSampleSize << " centers." << endl; ////
			for (unsigned i = 0; i < density_list.size() && result.size() < aSampleSize; i++) {
				unsigned id = density_list[i].second;
				vector<unsigned> neighborhood = ComputeApproximateNeighborhood(id);
				set<unsigned> neighborhood_set;
				neighborhood_set.insert(neighborhood.begin(), neighborhood.end());
				set<unsigned> intersection;
				set_intersection(active_neighborhood.begin(), active_neighborhood.end(), neighborhood_set.begin(), neighborhood_set.end(), inserter(intersection, intersection.begin()));
				if (i == 0 || intersection.size() <= aMaxIntersectionSize) { //if the intersection between the neighborhood of the current center and the union of all active neighborhoods is less than a defined constant (eg. 0) then accept the new center in the active set
					active_neighborhood.insert(neighborhood.begin(), neighborhood.end());
					result.push_back(id);
					progress_bar.Count();
				}
			}
		}
		//return result;

		//iterate the following steps:
		//swap one center
		//if it improves the overall quality measure keep it otherwise do not perform the update
		vector<unsigned> working_set(result);
		if (PARAM_OBJ.mMaxRefinement > 0 && PARAM_OBJ.mMaxIntersectionSize > 0) {
			cout << "Attempt to improve over initial result (max num iterations " << PARAM_OBJ.mMaxRefinement << ")" << endl;
			unsigned counter = 0;
			for (unsigned j = 0; j < PARAM_OBJ.mMaxRefinement; j++) {
				vector<unsigned> current_set(working_set);
				double original_quality = CenterSetQuality(current_set);
				unsigned rand_index_dest = randomUnsigned(selected_index_list.size());
				unsigned dest_id = selected_index_list[rand_index_dest];
				unsigned rand_index_src = randomUnsigned(current_set.size());
				unsigned src_id = current_set[rand_index_src];
				current_set[rand_index_src] = dest_id;
				double current_quality = CenterSetQuality(current_set);
				if (current_quality > original_quality) {
					working_set = current_set;
					counter++;
					cout << "Replaced center with id=" << src_id << " with id=" << dest_id << endl;
				}
			}
			if (counter > 0) {
				cout << "Successfully executed " << counter << " replacements" << endl;
			}
		}
		return working_set;
	}

	double CenterSetQuality(vector<unsigned> aSet) {
		//compute sum of approximate densities
		double density = 0;
		for (unsigned i = 0; i < aSet.size(); ++i) {
			density += ComputeApproximateDensity(aSet[i]);
		}
		//compute scaling factor as the ration of the size of the union of the neighb over the sum of the sizes of the neigh

		unsigned size_sum = 0;
		set<unsigned> union_set;
		for (unsigned i = 0; i < aSet.size(); ++i) {
			vector<unsigned> neighborhood = ComputeApproximateNeighborhood(aSet[i]);
			union_set.insert(neighborhood.begin(), neighborhood.end());
			size_sum += neighborhood.size();
		}
		double scaling_factor = (double) union_set.size() / (double) size_sum;
		double set_quality = density * scaling_factor;
		return set_quality;
	}

	double ComputeAverageSimilarityToSet(unsigned aID, set<unsigned>& aSet) {
		double avg_sim = 0;
		for (set<unsigned>::iterator it = aSet.begin(); it != aSet.end(); ++it) {
			unsigned id = (*it);
			avg_sim += Kernel(aID, id);
		}
		return avg_sim / aSet.size();
	}

	double ComputeAverageSimilarity(set<unsigned>& aSet) {
		double avg_sim = 0;
		for (set<unsigned>::iterator it = aSet.begin(); it != aSet.end(); ++it) {
			unsigned id = (*it);
			avg_sim += ComputeAverageSimilarityToSet(id, aSet);
		}
		return avg_sim / aSet.size();
	}

	vector<unsigned> ComputeNeighborhoodRanking(unsigned aID) {
		vector<pair<double, unsigned> > sim_list;
		for (unsigned i = 0; i < mDataset.size(); ++i) {
			if (i != aID) {
				double k = Kernel(aID, i);
				sim_list.push_back(make_pair(-k, i)); //note: use -k to sort in decreasing order
			}
		}
		//sort and take most similar <effective_size>
		sort(sim_list.begin(), sim_list.end());
		vector<unsigned> neighborhood;
		for (unsigned t = 0; t < sim_list.size(); ++t)
			neighborhood.push_back(sim_list[t].second);
		return neighborhood;
	}

	void OutputApproximateCluster() {
//TODO finish the approximate cluster
		//output cluster but replace all calls to true_similarity with shared neighb similarity over approximate neighb
	}

	void OutputCluster(ostream& out) {
		//initialize density cache
		mApproximateDensityMap.clear();
		mTrueDensityMap.clear();
		for (unsigned i = 0; i < mDataset.size(); ++i) {
			mApproximateDensityMap.push_back(-1);
			mTrueDensityMap.push_back(-1);
		}

		vector<unsigned> density_center_list;
		density_center_list = ComputeMinimallyOverlappingHighDensityCenterList(PARAM_OBJ.mSampleSize, PARAM_OBJ.mFractionCenterScan, PARAM_OBJ.mMaxIntersectionSize);
		cout << "Compute (approximate) neighborhood for selected " << density_center_list.size() << " cluster centers." << endl; ////
		{
			ProgressBar progress_bar;
			for (unsigned i = 0; i < density_center_list.size(); ++i) {
				unsigned id = density_center_list[i];
				vector<unsigned> neighborhood = ComputeTrueNeighborhood(id);
				progress_bar.Count();
				for (unsigned i = 0; i < neighborhood.size(); i++) {
					unsigned relative_id = neighborhood[i];
					unsigned absolute_id = mIdMap[relative_id];
					out << absolute_id << " ";
				}
				out << endl;
			}
		}
		if (PARAM_OBJ.mVerbose) OutputClusterVerbose(out, density_center_list);
	}

	void OutputClusterVerbose(ostream& out, vector<unsigned> aDensityCenterList) {
		vector<double> density_list;
		{
			ProgressBar progress_bar;
			cout << "Computing true density for approximated center list of " << aDensityCenterList.size() << " centers" << endl;
			for (unsigned i = 0; i < aDensityCenterList.size(); ++i) {
				unsigned id = aDensityCenterList[i];
				density_list.push_back(ComputeTrueDensity(id));
				progress_bar.Count();
			}
		}
		{
			VectorClass density_stats(density_list);
			cout << endl << "Centers density statistics: ";
			density_stats.OutputStatistics(cout);
			cout << endl;
		}
		{
			vector<double> true_density_list;
			{
				ProgressBar progress_bar;
				cout << "Computing true density for all instances " << mDataset.size() << endl;
				for (unsigned i = 0; i < mDataset.size(); i++) {
					true_density_list.push_back(ComputeTrueDensity(i));
					progress_bar.Count();
				}
			}
			VectorClass density_stats(true_density_list);
			cout << endl << "Global density statistics: ";
			density_stats.OutputStatistics(cout);
			cout << endl;
		}
	}

	void OutputAccuracy(ostream& out) {
		cout << "Compute neighbourhood accuracy" << endl; ////
		ProgressBar progress_bar;
		double cum = 0;
		unsigned effective_neighbourhood_size = min(PARAM_OBJ.mNumNearestNeighbors, (unsigned) mDataset.size());
		for (unsigned u = 0; u < mDataset.size(); ++u) {
			progress_bar.Count();
			vector<unsigned> approximate_neighborhood = ComputeApproximateNeighborhood(u);
			set<unsigned> approximate_neighborhood_set;
			approximate_neighborhood_set.insert(approximate_neighborhood.begin(), approximate_neighborhood.end());
			vector<unsigned> true_neighborhood = ComputeTrueNeighborhood(u);
			set<unsigned> true_neighborhood_set;
			true_neighborhood_set.insert(true_neighborhood.begin(), true_neighborhood.end());
			set<unsigned> intersection;
			set_intersection(approximate_neighborhood_set.begin(), approximate_neighborhood_set.end(), true_neighborhood_set.begin(), true_neighborhood_set.end(), inserter(intersection, intersection.begin()));
			double val = (double) intersection.size() / effective_neighbourhood_size;
			out << val << endl;
			cum += val;
		}
		cout << endl << "Accuracy: " << cum / mDataset.size() << endl;
	}

	void OutputApproximateKNN(ostream& out) {
		vector<int> id_list;
		if (mGreyList.size() > 0) {
			for (unsigned i = 0; i < mGreyList.size(); ++i)
				id_list.push_back(mGreyList[i] - 1);
		} else {
			for (unsigned i = 0; i < mDataset.size(); ++i)
				id_list.push_back(i);
		}
		cout << "Compute approximate nearest neighbours for " << id_list.size() << " elements." << endl; ////
		ProgressBar progress_bar;
		for (unsigned u = 0; u < id_list.size(); ++u) {
			vector<unsigned> approximate_neighborhood = ComputeApproximateNeighborhood(id_list[u]);
			for (unsigned t = 0; t < approximate_neighborhood.size(); ++t) {
				out << approximate_neighborhood[t] + 1 << " "; //NOTE: numbering starts from 1
			}
			out << endl;
			progress_bar.Count();
		}
	}

	void OutputTrueKNN(ostream& out) {
		vector<int> id_list;
		if (mGreyList.size() > 0) {
			for (unsigned i = 0; i < mGreyList.size(); ++i)
				id_list.push_back(mGreyList[i] - 1);
		} else {
			for (unsigned i = 0; i < mDataset.size(); ++i)
				id_list.push_back(i);
		}
		unsigned effective_neighbourhood_size = min(PARAM_OBJ.mNumNearestNeighbors, (unsigned) mDataset.size());
		cout << "Compute true " << effective_neighbourhood_size << "-nearest neighbours for " << id_list.size() << " elements." << endl; ////
		ProgressBar progress_bar;
		for (unsigned i = 0; i < id_list.size(); ++i) {
			unsigned u = id_list[i];
			vector<pair<double, unsigned> > sim_list;
			//compute kernel pairs between aID and all elements
			for (unsigned v = 0; v < mDataset.size(); v++) {
				double k_uv = Kernel(u, v);
				sim_list.push_back(make_pair(-k_uv, v)); //note: use -k to sort in decreasing order
			}
			//sort and take truly most similar
			sort(sim_list.begin(), sim_list.end());
			for (unsigned k = 0; k < effective_neighbourhood_size; ++k) {
				out << sim_list[k].second + 1 << " "; //NOTE: numbering starts from 1
			}
			out << endl;
			progress_bar.Count();
		}
	}

	void OutputApproximateKNNPrediction(ostream& out, string& aSparseASCIITestInputFileName, string& aTrainTargetInputFileName) {
		cout << "Compute approximate nearest neighbour prediction of test instances in " << aSparseASCIITestInputFileName << "." << endl; ////
		ProgressBar progress_bar;

		//read test instances
		vector<SVector> test_dataset;
		InputSparse(aSparseASCIITestInputFileName, "ascii", test_dataset);

		//read train targets
		vector<string> train_target_list;
		InputStringList(aTrainTargetInputFileName, train_target_list);

		cout << "Computing k-NN predictions." << endl;
		//for each test instance
		for (unsigned u = 0; u < test_dataset.size(); ++u) {
			//extract signature
			vector<unsigned> hash_signature = ComputeHashSignature(test_dataset[u], u);
			//extract knn
			vector<unsigned> approximate_neighborhood = ComputeApproximateNeighborhood(hash_signature);
			string prediction = KNNPredict(approximate_neighborhood, train_target_list);
			out << prediction << endl;
			progress_bar.Count();
		}
	}

	void OutputTrueKNNPrediction(ostream& out, string& aSparseASCIITestInputFileName, string& aTrainTargetInputFileName) {
		cout << "Compute true nearest neighbour prediction of test instances." << endl; ////
		ProgressBar progress_bar;

		//read test instances
		vector<SVector> test_dataset;
		InputSparse(aSparseASCIITestInputFileName, "ascii", test_dataset);

		//read train targets
		vector<string> train_target_list;
		InputStringList(aTrainTargetInputFileName, train_target_list);

		cout << "Computing k-NN predictions on " << test_dataset.size() << " instances." << endl;
		//for each test instance
		for (unsigned u = 0; u < test_dataset.size(); ++u) {
			//extract knn
			vector<unsigned> neighborhood = ComputeTrueNeighborhood(test_dataset[u]);
			//compute majority class
			string prediction = KNNPredict(neighborhood, train_target_list);
			out << prediction << endl;
			progress_bar.Count();
		}
	}

	string KNNPredict(const vector<unsigned>& aNeighborhood, const vector<string>& aTargetList) const {
		//compute histogram of targets in neighborhood
		map<string, unsigned> histogram;
		for (unsigned i = 0; i < aNeighborhood.size(); ++i) {
			unsigned nn_id = aNeighborhood[i];
			assert(nn_id<aTargetList.size());
			string predicted_target = aTargetList[nn_id];
			if (histogram.count(predicted_target) == 0) histogram[predicted_target] = 1;
			else histogram[predicted_target]++;
		}
		//compute majority vote for target
		string max_target = aTargetList[0]; //initialization with one arbitrary target
		unsigned max_val = 0;
		for (map<string, unsigned>::const_iterator it = histogram.begin(); it != histogram.end(); ++it) {
			string target = it->first;
			unsigned vote = it->second;
			if (max_val < vote) {
				max_val = vote;
				max_target = target;
			}
		}
		return max_target;
	}

	void OutputKernel(ostream& out) {
		cout << "Compute kernel matrix." << endl; ////
		ProgressBar progress_bar;
		for (unsigned i = 0; i < mDataset.size(); i++) {
			for (unsigned j = 0; j < mDataset.size(); j++)
				out << Kernel(i, j) << " ";
			out << endl;
			progress_bar.Count();
		}
	}

	void Output(ostream& out) {
		for (unsigned i = 0; i < mDataset.size(); i++)
			out << mDataset[i];
	}

	void OutputFeatureMap(string aFileName) const {
		string ofname = aFileName + ".feature_map";
		ofstream of(ofname.c_str());
		pmFeatureGenerator->OutputFeatureMap(of);
	}

	double Kernel(unsigned aI, unsigned aJ) {
		unsigned i = min(aI, aJ);
		unsigned j = max(aI, aJ);
		if (PARAM_OBJ.mCaching) {
			pair<unsigned, unsigned> key = make_pair(i, j);
			if (mKernelMap.count(key) == 0) {
				double value = Similarity(i, j);
				mKernelMap[key] = value;
			}
			return mKernelMap[key];
		} else return Similarity(i, j);
	}

	double Similarity(unsigned aI, unsigned aJ) {
		return dot(mDataset[aI], mDataset[aJ]);
	}
	/**
	 Computes the fraction of neighbors that are common between instance I and J
	 */
	double ComputeSharedNeighborhoodSimilarity(unsigned aI, unsigned aJ) {
		//TODO: cache this result
		vector<unsigned> neighborhood_i = ComputeApproximateNeighborhood(aI);
		vector<unsigned> neighborhood_j = ComputeApproximateNeighborhood(aJ);
		set<unsigned> neighborhood_i_set;
		neighborhood_i_set.insert(neighborhood_i.begin(), neighborhood_i.end());
		set<unsigned> neighborhood_j_set;
		neighborhood_j_set.insert(neighborhood_j.begin(), neighborhood_j.end());
		set<unsigned> intersection;
		set_intersection(neighborhood_i_set.begin(), neighborhood_i_set.end(), neighborhood_j_set.begin(), neighborhood_j_set.end(), inserter(intersection, intersection.begin()));
		double shared_neighbourhood_value = (double) intersection.size() / sqrt((double) neighborhood_i_set.size() * (double) neighborhood_j_set.size());
		return shared_neighbourhood_value;
	}
}
;

//---------------------------------------------------------------------------------
int main(int argc, char** argv) {
	TimerClass T;
	LOGF << "--------------------------------------------------------------------------------" << endl;
	LOGF << CREDITS << " \nLast build:" << DATE << endl;
	cout << CREDITS << " \nLast build:" << DATE << endl;
	time_t rawtime_start;
	time(&rawtime_start);
	LOGF << "Start logging: " << asctime(localtime(&rawtime_start)) << endl;
	LOGF << "Command line: ";
	for (int i = 0; i < argc; i++)
		LOGF << stream_cast<string>(argv[i]) << " ";
	LOGF << endl;
	try {
		PARAM_OBJ.Init(argc, argv);
		srand(PARAM_OBJ.mRandSeed);
		string mode = "";
		//factory
		NSPDK_FeatureGenerator fg("nspdk");

		NSPDK_FeatureGenerator* pfg;
		if (PARAM_OBJ.mType == "nspdk") pfg = &fg;
		else throw range_error("Unknown feature generator type:" + PARAM_OBJ.mType);

		if (PARAM_OBJ.mOutputFeatureMap) PARAM_OBJ.mDebug += 1; //if the output of the feature map is required then the Debug level has to be at least 1

		pfg->set_flag("radius", stream_cast<string>(PARAM_OBJ.mRadiusMax));
		pfg->set_flag("distance", stream_cast<string>(PARAM_OBJ.mDistanceMax));
		pfg->set_flag("match_type", stream_cast<string>(PARAM_OBJ.mMatchingType));
		pfg->set_flag("hash_bit_size", stream_cast<string>(PARAM_OBJ.mFeatureBitSize));
		pfg->set_flag("hash_bit_mask", stream_cast<string>((2 << PARAM_OBJ.mFeatureBitSize) - 1));
		pfg->set_flag("verbosity", stream_cast<string>(PARAM_OBJ.mDebug));
		if (PARAM_OBJ.mMinKernel) pfg->set_flag("min_kernel", "true");
		if (!PARAM_OBJ.mNormalization) pfg->set_flag("normalization", "false");
		pfg->OutputParameters(cout); ////////////////////////////////////////////////////////////////

		string ofname;
		//main process
		NSPDKClass C(pfg);

		if (PARAM_OBJ.mOutputApproximateCluster) C.InitBinDataStructure();

		//read data
		C.Input();

		if (PARAM_OBJ.mOutputFeatures) {
			if (PARAM_OBJ.mOutputFeatureMap) C.OutputFeatureMap(PARAM_OBJ.mGspanInputFileName);
		}

		if (PARAM_OBJ.mOutputApproximateCluster) {
			C.OutputApproximateCluster();
		}

		if (PARAM_OBJ.mOutputCluster) {
			C.ComputeBinDataStructure();
			if (PARAM_OBJ.mVerbose) C.OutputBinDataStructureStatistics();
			ofname = PARAM_OBJ.mGspanInputFileName + PARAM_OBJ.mSparseASCIIInputFileName + PARAM_OBJ.mSparseBinaryInputFileName + ".fast_cluster";
			ofstream ofs_fc(ofname.c_str());
			C.OutputCluster(ofs_fc);
		}

		if (PARAM_OBJ.mOutputAccuracy) {
			C.ComputeBinDataStructure();
			ofname = PARAM_OBJ.mGspanInputFileName + PARAM_OBJ.mSparseASCIIInputFileName + PARAM_OBJ.mSparseBinaryInputFileName + ".nn_accuracy";
			ofstream ofs_fc(ofname.c_str());
			C.OutputAccuracy(ofs_fc);
		}

		if (PARAM_OBJ.mOutputApproximateKNN) {
			C.ComputeBinDataStructure();
			ofname = PARAM_OBJ.mGspanInputFileName + PARAM_OBJ.mSparseASCIIInputFileName + PARAM_OBJ.mSparseBinaryInputFileName + ".approx_knn";
			ofstream ofs_aknn(ofname.c_str());
			C.OutputApproximateKNN(ofs_aknn);
			cout << endl << "Results written in file <" << ofname << ">" << endl;
		}

		if (PARAM_OBJ.mOutputTrueKNN) {
			ofname = PARAM_OBJ.mGspanInputFileName + PARAM_OBJ.mSparseASCIIInputFileName + PARAM_OBJ.mSparseBinaryInputFileName + ".knn";
			ofstream ofs_fknn(ofname.c_str());
			C.OutputTrueKNN(ofs_fknn);
		}

		if (PARAM_OBJ.mOutputKernel) {
			ofname = PARAM_OBJ.mGspanInputFileName + PARAM_OBJ.mSparseASCIIInputFileName + PARAM_OBJ.mSparseBinaryInputFileName + ".kernel";
			ofstream ofs_fk(ofname.c_str());
			C.OutputKernel(ofs_fk);
		}

		if (PARAM_OBJ.mOutputApproximateKNNPrediction) {
			C.ComputeBinDataStructure();
			ofname = PARAM_OBJ.mGspanInputFileName + PARAM_OBJ.mSparseASCIIInputFileName + PARAM_OBJ.mSparseBinaryInputFileName + ".approx_knn_prediction";
			ofstream ofs_knnp(ofname.c_str());
			C.OutputApproximateKNNPrediction(ofs_knnp, PARAM_OBJ.mSparseASCIITestInputFileName, PARAM_OBJ.mTrainTargetInputFileName);
		}

		if (PARAM_OBJ.mOutputTrueKNNPrediction) {
			ofname = PARAM_OBJ.mGspanInputFileName + PARAM_OBJ.mSparseASCIIInputFileName + PARAM_OBJ.mSparseBinaryInputFileName + ".knn_prediction";
			ofstream ofs_knnp(ofname.c_str());
			C.OutputTrueKNNPrediction(ofs_knnp, PARAM_OBJ.mSparseASCIITestInputFileName, PARAM_OBJ.mTrainTargetInputFileName);
		}

	} catch (exception& e) {
		cerr << e.what();
		LOGF << e.what() << endl;
	}
	time_t rawtime_end;
	time(&rawtime_end);
	LOGF << "End logging: " << asctime(localtime(&rawtime_end)) << endl;
	return 0;
}
