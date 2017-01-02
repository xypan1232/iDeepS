#include "Utility.h"
#include "BaseGraphClass.h"
#include "GraphClass.h"
#include "NSPDK_FeatureGenerator.h"
#include "vectors.h"
#include "gzstream.h"
#include "OpenBabelConverter.h"
#include "DDK_FeatureGenerator.h"

#include <omp.h>

using namespace std;

const string PROG_CREDIT = "EDeN (Explicit Decomposition with Neighborhoods) Vers. 0.4.2 (15 June 2013)\nAuthor: Fabrizio Costa costa@informatik.uni-freiburg.de";
const string CITATIONS =
		"For the Decompositional DAG graph kernel see G. Da San Martino, N. Navarin and A. Sperduti , ''A tree-based kernel for graphs'', Proceedings of the Twelfth SIAM International Conference on Data Mining, Anaheim, California, April 26 - 28, 2012, p. 975-986.\n"
		"For the NSPDK graph kernels see Fabrizio Costa, Kurt De Grave, ''Fast Neighborhood Subgraph Pairwise Distance Kernel'', Proceedings of the 27th International Conference on Machine Learning (ICML-2010), Haifa, Israel, 2010..\n"
		"The code for Stochastic Gradient Descent SVM is adapted from http://leon.bottou.org/projects/sgd. Léon Bottou and Yann LeCun, ''Large Scale Online Learning'', Advances in Neural Information Processing Systems 16, Edited by Sebastian Thrun, Lawrence Saul and Bernhard Schölkopf, MIT Press, Cambridge, MA, 2004.\n"
		"The embedding method is adapted from L.Chen, A.Buja ''Local Multidimensional Scaling for Nonlinear Dimension Reduction, Graph Drawing, and Proximity Analysis'', Journal of the American Statistical Association, 2009.";

FlagsService& The_FlagsService = FlagsService::get_instance();

enum ActionType {
	NULL_ACTION, TRAIN, TEST, PARAMETERS_OPTIMIZATION, CROSS_VALIDATION, BIAS_VARIANCE_DECOMPOSITION, LEARNING_CURVE, TEST_PART, FEATURE, FEATURE_PART, FEATURE_SCALED, MATRIX, EMBED, TARGET_ALIGNMENT, CLUSTER, NEAREST_NEIGHBOR, MIN_HASH, SEMI_SUPERVISED
};

enum InputFileType {
	SPARSE_VECTOR, GRAPH, MOLECULAR_GRAPH, SEQUENCE
};

const unsigned MAXUNSIGNED = 2 << 30;
const unsigned BUFFER_SIZE = 500;

//------------------------------------------------------------------------------------------------------------------------
// Available losses
#define HINGELOSS 1
#define SMOOTHHINGELOSS 2
#define SQUAREDHINGELOSS 3
#define LOGLOSS 10
#define LOGLOSSMARGIN 11

// Select loss: NOTE: the selection done in the makefile
//#define LOSS LOGLOSS
//#define LOSS HINGELOSS
//#define LOSS SMOOTHHINGELOSS

// Zero when no bias
// One when bias term
#define BIAS 1

inline
double loss(double z) {
#if LOSS == LOGLOSS
	if (z > 18)
	return exp(-z);
	if (z < -18)
	return -z;
	return log(1+exp(-z));
#elif LOSS == LOGLOSSMARGIN
	if (z > 18)
	return exp(1-z);
	if (z < -18)
	return 1-z;
	return log(1+exp(1-z));
#elif LOSS == SMOOTHHINGELOSS
	if (z < 0)
	return 0.5 - z;
	if (z < 1)
	return 0.5 * (1-z) * (1-z);
	return 0;
#elif LOSS == SQUAREDHINGELOSS
	if (z < 1)
	return 0.5 * (1 - z) * (1 - z);
	return 0;
#elif LOSS == HINGELOSS
	if (z < 1)
		return 1 - z;
	return 0;
#else
	return 0;
#endif
}

inline
double dloss(double z) {
#if LOSS == LOGLOSS
	if (z > 18)
	return exp(-z);
	if (z < -18)
	return 1;
	return 1 / (exp(z) + 1);
#elif LOSS == LOGLOSSMARGIN
	if (z > 18)
	return exp(1-z);
	if (z < -18)
	return 1;
	return 1 / (exp(z-1) + 1);
#elif LOSS == SMOOTHHINGELOSS
	if (z < 0)
	return 1;
	if (z < 1)
	return 1-z;
	return 0;
#elif LOSS == SQUAREDHINGELOSS
	if (z < 1)
	return (1 - z);
	return 0;
#else
	if (z < 1)
		return 1;
	return 0;
#endif
}

//------------------------------------------------------------------------------------------------------------------------
class Parameters {
public:
	map<string, ParameterType> mOptionList;
	map<ActionType, vector<ParameterType*> > mOptionGrouping;

	string mAction;
	ActionType mActionCode;
	string mInputDataFileName;
	string mTargetFileName;
	string mRowIndexFileName;
	string mColIndexFileName;
	string mModelFileName;
	string mFileType;
	bool mBinaryFormat;
	string mOpenBabelFormat;
	InputFileType mFileTypeCode;
	bool mKernelNoNormalization;
	bool mMinKernel;
	unsigned mRadius;
	unsigned mDistance;
	unsigned mVertexDegreeThreshold;
	double mLambda;
	unsigned mEpochs;
	unsigned mHashBitSize;
	unsigned mCrossValidationNumFolds;
	unsigned mNumPoints;
	unsigned mRandomSeed;
	string mKernelType;
	string mGraphType;
	unsigned mSemiSupervisedNumIterations;
	double mSemiSupervisedThreshold;
	bool mSemiSupervisedInduceOnlyPositive;
	bool mSemiSupervisedInduceOnlyNegative;
	string mSuffix;
	unsigned mSequenceDegree;
	unsigned mLMDSNumRandomRestarts;
	unsigned mLMDSNumIterations;
	unsigned mLMDSDimensionality;
	double mLMDSIterationEpsilon;
	unsigned mLMDSNeighborhoodSize;
	unsigned mLMDSNonNeighborhoodSize;
	unsigned mLMDSNeighborhoodSizeRange;
	double mLMDSTau;
	unsigned mLMDSTauExponentRange;
	bool mVerbose;
	bool mMinimalOutput;
	bool mSequenceToken;
	bool mSequenceMultiLine;
	bool mSequencePairwiseInteraction;
	unsigned mSparsificationNumIterations;
	unsigned mTopologicalRegularizationNumNeighbors;
	double mTopologicalRegularizationRate;
	unsigned mNumLineSearchIterations;

	unsigned mNumHashFunctions;
	unsigned mNumRepeatsHashFunction;
	double mMaxSizeBin;
	double mEccessNeighbourSizeFactor;
	unsigned mSampleSize;
	unsigned mNumNearestNeighbors;
	bool mSharedNeighborhood;
	double mFractionCenterScan;
	unsigned mMaxIntersectionSize;
	double mClusterThreshold;
	string mClusterType;
	bool mNoNeighborhoodCache;
	bool mNoMinHashCache;
	bool mForceApproximate;

	double mSemiSupervisedAlpha;

	//DD kernel family additional parameters
	double mTreeLambda;
	double mRadiusTwo;

public:
	Parameters() {
		SetupOptions();
	}

	void SetupOptions() {
		{
			ParameterType param;
			param.mShortSwitch = "h";
			param.mLongSwitch = "help";
			param.mShortDescription = "Prints compact help.";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
		}
		{
			ParameterType param;
			param.mShortSwitch = "H";
			param.mLongSwitch = "Help";
			param.mShortDescription = "Prints extended help.";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
		}
		{
			ParameterType param;
			param.mShortSwitch = "a";
			param.mLongSwitch = "action";
			param.mShortDescription = "";
			param.mTypeCode = LIST;
			param.mValue = "";
			param.mCloseValuesList.push_back("TRAIN");
			param.mCloseValuesList.push_back("TEST");
			param.mCloseValuesList.push_back("TEST_PART");
			param.mCloseValuesList.push_back("CROSS_VALIDATION");
			param.mCloseValuesList.push_back("BIAS_VARIANCE_DECOMPOSITION");
			param.mCloseValuesList.push_back("PARAMETERS_OPTIMIZATION");
			param.mCloseValuesList.push_back("LEARNING_CURVE");
			param.mCloseValuesList.push_back("FEATURE");
			param.mCloseValuesList.push_back("FEATURE_PART");
			param.mCloseValuesList.push_back("FEATURE_SCALED");
			param.mCloseValuesList.push_back("MATRIX");
			param.mCloseValuesList.push_back("EMBED");
			param.mCloseValuesList.push_back("TARGET_ALIGNMENT");
			param.mCloseValuesList.push_back("CLUSTER");
			param.mCloseValuesList.push_back("NEAREST_NEIGHBOR");
			param.mCloseValuesList.push_back("MIN_HASH");
			param.mCloseValuesList.push_back("SEMI_SUPERVISED");

			mOptionList.insert(make_pair(param.mLongSwitch, param));

			mOptionGrouping.insert(make_pair(TRAIN, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(TEST, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(TEST_PART, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(CROSS_VALIDATION, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(BIAS_VARIANCE_DECOMPOSITION, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(PARAMETERS_OPTIMIZATION, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(LEARNING_CURVE, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(FEATURE, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(FEATURE_PART, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(FEATURE_SCALED, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(MATRIX, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(EMBED, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(TARGET_ALIGNMENT, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(CLUSTER, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(NEAREST_NEIGHBOR, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(MIN_HASH, vector<ParameterType*>()));
			mOptionGrouping.insert(make_pair(SEMI_SUPERVISED, vector<ParameterType*>()));
		}
		{
			ParameterType param;
			param.mShortSwitch = "i";
			param.mLongSwitch = "input_data_file_name";
			param.mShortDescription = "";
			param.mTypeCode = STRING;
			param.mValue = "";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "t";
			param.mLongSwitch = "target_file_name";
			param.mShortDescription = "";
			param.mTypeCode = STRING;
			param.mValue = "";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "A";
			param.mLongSwitch = "row_index_file_name";
			param.mShortDescription = "In MATRIX and in NEAREST_NEIGHBOR.";
			param.mTypeCode = STRING;
			param.mValue = "";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "B";
			param.mLongSwitch = "col_index_file_name";
			param.mShortDescription = "In MATRIX and in NEAREST_NEIGHBOR.";
			param.mTypeCode = STRING;
			param.mValue = "";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "f";
			param.mLongSwitch = "file_type";
			param.mShortDescription = "";
			param.mTypeCode = LIST;
			param.mValue = "GRAPH";
			param.mCloseValuesList.push_back("SPARSE_VECTOR");
			param.mCloseValuesList.push_back("GRAPH");
#ifdef USEOBABEL
			param.mCloseValuesList.push_back("MOLECULAR_GRAPH");
#endif
			param.mCloseValuesList.push_back("SEQUENCE");

			mOptionList.insert(make_pair(param.mLongSwitch, param));

			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}

		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "binary_file_type";
			param.mShortDescription = "";
			param.mTypeCode = FLAG;
			param.mValue = "0";

			mOptionList.insert(make_pair(param.mLongSwitch, param));

			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "k";
			param.mLongSwitch = "kernel_type";
			param.mShortDescription = "";
			param.mTypeCode = LIST;
			param.mValue = "NSPDK";
			param.mCloseValuesList.push_back("NSPDK");
			param.mCloseValuesList.push_back("WDK");
			param.mCloseValuesList.push_back("PBK");
			param.mCloseValuesList.push_back("USPK");
			param.mCloseValuesList.push_back("DDK");
			param.mCloseValuesList.push_back("NSDDK");
			param.mCloseValuesList.push_back("ANSDDK");


			param.mCloseValuesList.push_back("SK");

			mOptionList.insert(make_pair(param.mLongSwitch, param));

			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "g";
			param.mLongSwitch = "graph_type";
			param.mShortDescription = "";
			param.mTypeCode = LIST;
			param.mValue = "UNDIRECTED";
			param.mCloseValuesList.push_back("DIRECTED");
			param.mCloseValuesList.push_back("UNDIRECTED");
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "no_normalization";
			param.mShortDescription = "Kernel parameter.";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "min_kernel";
			param.mShortDescription = "Kernel parameter.";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
#ifdef USEOBABEL
		{
			ParameterType param;
			param.mShortSwitch = "o";
			param.mLongSwitch = "open_babel_file_format";
			param.mShortDescription = "MOLECULAR_GRAPH parameter.";
			param.mTypeCode = STRING;
			param.mValue = "sdf";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
#endif
		{
			ParameterType param;
			param.mShortSwitch = "m";
			param.mLongSwitch = "model_file_name";
			param.mShortDescription = "";
			param.mTypeCode = STRING;
			param.mValue = "model";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "s";
			param.mLongSwitch = "suffix";
			param.mShortDescription = "Suffix string for all output files.";
			param.mTypeCode = STRING;
			param.mValue = "";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "b";
			param.mLongSwitch = "hash_bit_size";
			param.mShortDescription = "Kernel parameter.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "15";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "r";
			param.mLongSwitch = "radius";
			param.mShortDescription = "Kernel parameter.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "2";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "d";
			param.mLongSwitch = "distance";
			param.mShortDescription = "Kernel parameter.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "5";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "v";
			param.mLongSwitch = "vertex_degree_threshold";
			param.mShortDescription = "Kernel parameter.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "7";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "M";
			param.mLongSwitch = "sequence_degree";
			param.mShortDescription = "SEQUENCE data type.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "1";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "sequence_token";
			param.mShortDescription = "Labels are strings separated by spaces.";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "sequence_multi_line";
			param.mShortDescription = "The annotation is encoded on subsequent rows.";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "sequence_pairwise_interaction";
			param.mShortDescription = "Abstraction vertices relating all pairs of disjointed sequences are inserted.";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "l";
			param.mLongSwitch = "lambda";
			param.mShortDescription = "Stochastic gradient descend algorithm.";
			param.mTypeCode = REAL;
			param.mValue = "1e-4";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "e";
			param.mLongSwitch = "epochs";
			param.mShortDescription = "Stochastic gradient descend algorithm.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "10";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "O";
			param.mLongSwitch = "sparsification_num_iterations";
			param.mShortDescription = "In training.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "C";
			param.mLongSwitch = "topological_regularization_num_neighbors";
			param.mShortDescription = "In training.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "L";
			param.mLongSwitch = "topological_regularization_decay_rate";
			param.mShortDescription = "In training.";
			param.mTypeCode = REAL;
			param.mValue = "0.01";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "c";
			param.mLongSwitch = "num_cross_validation_folds";
			param.mShortDescription = "";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "10";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "p";
			param.mLongSwitch = "num_evaluation_points";
			param.mShortDescription = "";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "10";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "Y";
			param.mLongSwitch = "num_line_search_iterations";
			param.mShortDescription = "";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "3";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
		}
		{
			ParameterType param;
			param.mShortSwitch = "S";
			param.mLongSwitch = "num_iterations";
			param.mShortDescription = "In semi-supervised setting.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "3";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "T";
			param.mLongSwitch = "threshold";
			param.mShortDescription = "In semi-supervised setting. Only the top and low quantile will be used as positives and negative instances. A threshold of 1 means that all unsupervised instaces are used in the next phase.";
			param.mTypeCode = REAL;
			param.mValue = "1";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "only_positive";
			param.mShortDescription = "In semi-supervised setting. Induce only positive class instances.";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "only_negative";
			param.mShortDescription = "In semi-supervised setting. Induce only negative class instances.";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "N";
			param.mLongSwitch = "num_of_random_restarts";
			param.mShortDescription = "In EMBED.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "1";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "I";
			param.mLongSwitch = "num_of_iterations";
			param.mShortDescription = "In EMBED.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "5000";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "D";
			param.mLongSwitch = "dimensionality";
			param.mShortDescription = "In EMBED.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "2";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "E";
			param.mLongSwitch = "epsilon";
			param.mShortDescription = "In EMBED.";
			param.mTypeCode = REAL;
			param.mValue = "0.01";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "K";
			param.mLongSwitch = "neighborhood_size";
			param.mShortDescription = "In EMBED.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "10";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "n";
			param.mLongSwitch = "non_neighborhood_size";
			param.mShortDescription = "In EMBED.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "100";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "G";
			param.mLongSwitch = "neighborhood_size_range";
			param.mShortDescription = "In EMBED.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "U";
			param.mLongSwitch = "tau";
			param.mShortDescription = "In EMBED.";
			param.mTypeCode = REAL;
			param.mValue = "0.0005";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "u";
			param.mLongSwitch = "tau_exponent_range";
			param.mShortDescription = "In EMBED.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "1";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "V";
			param.mLongSwitch = "verbose";
			param.mShortDescription = "Outputs the graphs and the feature encodings.";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "W";
			param.mLongSwitch = "minimal_output";
			param.mShortDescription = "Output only results.";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "cluster_type";
			param.mShortDescription = "";
			param.mTypeCode = LIST;
			param.mValue = "K_QUICK_SHIFT";
			param.mCloseValuesList.push_back("K_QUICK_SHIFT");
			param.mCloseValuesList.push_back("DENSE_CENTERS");
			param.mCloseValuesList.push_back("APPROXIMATION_ACCURACY");
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "force_approximate";
			param.mShortDescription = "In NEAREST_NEIGHBOUR";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "F";
			param.mLongSwitch = "num_hash_functions";
			param.mShortDescription = "In CLUSTER";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "400";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "num_repeat_hash_functions";
			param.mShortDescription = "In approximate neighborhood";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "10";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "z";
			param.mLongSwitch = "max_size_bin";
			param.mShortDescription =
					"In CLUSTER.  Expressed as the maximum fraction of the datset size. When a bin contains references to more instances than this quantity, then the bin is erased. The ratio is that this featrue is common to too many instances and it is therefore not informative. Morover the runtimes become non sub-linear if a significant fraction of the dataset size has to be checked.";
			param.mTypeCode = REAL;
			param.mValue = "0.1";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "eccess_neighbour_size_factor";
			param.mShortDescription = "In CLUSTER. Expressed as a multiplicative factor w.r.t. the neighbourhood size required. It means that the approximate neighborhood query stops at the X most frequent instances, where X= eccess_neighbour_size_factor * neighbourhood size.";
			param.mTypeCode = REAL;
			param.mValue = "5.0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "Q";
			param.mLongSwitch = "sample_size";
			param.mShortDescription = "In CLUSTER";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "20";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "x";
			param.mLongSwitch = "num_nearest_neighbours";
			param.mShortDescription = "In CLUSTER and in NEAREST_NEIGHBOR";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "10";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "X";
			param.mLongSwitch = "shared_neighborhood";
			param.mShortDescription = "In CLUSTER";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "fraction_center_scan";
			param.mShortDescription = "In CLUSTER";
			param.mTypeCode = REAL;
			param.mValue = "0.5";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "w";
			param.mLongSwitch = "max_intersection_size";
			param.mShortDescription = "In CLUSTER";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "3";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "cluster_threshold";
			param.mShortDescription = "Maximum factor for parent similarity w.r.t. average neighbors similarity (normalized). A smaller value (e.g. 0.5) implies a stricter similarity criterion. In CLUSTER";
			param.mTypeCode = REAL;
			param.mValue = "1";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "semi_supervised_alpha";
			param.mShortDescription = "In SEMI_SUPERVISED";
			param.mTypeCode = REAL;
			param.mValue = ".99";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "no_neighborhood_cache";
			param.mShortDescription = "In CLUSTER";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
		}
		{
			ParameterType param;
			param.mShortSwitch = "";
			param.mLongSwitch = "no_minhash_cache";
			param.mShortDescription = "In CLUSTER";
			param.mTypeCode = FLAG;
			param.mValue = "0";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
		}
		//DD kernel family additional parameters
		{
			ParameterType param;
			param.mLongSwitch = "tree_lambda";
			param.mShortDescription = "Kernel parameter.";
			param.mTypeCode = REAL;
			param.mValue = "1.2";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		{
			ParameterType param;
			param.mLongSwitch = "radius_two";
			param.mShortDescription = "Kernel parameter.";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "2";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TEST_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_PART];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[FEATURE_SCALED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[MATRIX];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[TARGET_ALIGNMENT];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CLUSTER];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[NEAREST_NEIGHBOR];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
		//--------------------------------------
		{
			ParameterType param;
			param.mShortSwitch = "R";
			param.mLongSwitch = "random_seed";
			param.mShortDescription = "";
			param.mTypeCode = POSITIVE_INTEGER;
			param.mValue = "1";
			mOptionList.insert(make_pair(param.mLongSwitch, param));
			{
				vector<ParameterType*>& vec = mOptionGrouping[TRAIN];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[CROSS_VALIDATION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[BIAS_VARIANCE_DECOMPOSITION];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[LEARNING_CURVE];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
			{
				vector<ParameterType*>& vec = mOptionGrouping[EMBED];
				ParameterType& p = mOptionList[param.mLongSwitch];
				vec.push_back(&p);
			}
		}
	}

	void Usage(string aCommandName, string aCompactOrExtended) {
		cout << SEP << endl << PROG_CREDIT << endl << SEP << endl;
		cout << "PROGRAM: " << aCommandName << endl;
		if (mOptionList["action"].mIsSet == false) {
			cout << "OPTIONS:" << endl;

			if (aCompactOrExtended == "EXTENDED")
				mOptionList["action"].OutputExtended(cout);
			else
				mOptionList["action"].OutputCompact(cout);
		} else {
			cout << "ACTION: " << mAction << endl;
			cout << "OPTIONS:" << endl;

			vector<ParameterType*> & vec = mOptionGrouping[mActionCode];
			for (unsigned i = 0; i < vec.size(); ++i) {
				assert(vec[i]);
				ParameterType& param = (*vec[i]);
				if (aCompactOrExtended == "EXTENDED")
					param.OutputExtended(cout);
				else
					param.OutputCompact(cout);
			}
		}
		cout << SEP << endl << "REFERENCES:" << endl << CITATIONS << endl << SEP << endl;
		exit(0);
	}

	/*
	 void UsageOld(string aCommandName, string aCompactOrExtended) {
	 cout << "PROGRAM: " << aCommandName << endl;
	 cout << "OPTIONS:" << endl;
	 for (unsigned i = 0; i < mOptionList.size(); ++i)
	 if (aCompactOrExtended == "EXTENDED")
	 mOptionList[i].OutputExtended(cout);
	 else
	 mOptionList[i].OutputCompact(cout);
	 cout << SEP << endl << CITATIONS << endl << SEP << endl;
	 exit(0);
	 }
	 */

	void Init(int argc, const char** argv) {
		if (argc == 1) {
			cout << "Use -h for compact help and -H for extended help." << endl;
			exit(1);
		}

		//convert argc in an option string vector
		vector<string> options;
		for (int i = 1; i < argc; i++)
			options.push_back(argv[i]);

		//parse the option string vector
		for (map<string, ParameterType>::iterator it = mOptionList.begin(); it != mOptionList.end(); ++it) {
			ParameterType& param = it->second;
			param.Parse(options);
		}

		//set the boolean parameters to a default value of false
		mKernelNoNormalization = false;
		mMinKernel = false;
		mSemiSupervisedInduceOnlyPositive = false;
		mSemiSupervisedInduceOnlyNegative = false;
		mVerbose = false;
		mSequenceToken = false;
		mSequenceMultiLine = false;
		mSequencePairwiseInteraction = false;
		mMinimalOutput = false;
		mSharedNeighborhood = false;
		mBinaryFormat = false;
		mForceApproximate = false;
		mNoNeighborhoodCache = false;
		mNoMinHashCache = false;

		//set the data members of Parameters according to user choice
		for (map<string, ParameterType>::iterator it = mOptionList.begin(); it != mOptionList.end(); ++it) {
			ParameterType& param = it->second;
			if (param.mIsSet) {
				if (param.mShortSwitch == "i")
					mInputDataFileName = param.mValue;
				if (param.mShortSwitch == "t")
					mTargetFileName = param.mValue;
				if (param.mLongSwitch == "no_normalization")
					mKernelNoNormalization = true;
				if (param.mLongSwitch == "min_kernel")
					mMinKernel = true;
				if (param.mLongSwitch == "only_positive")
					mSemiSupervisedInduceOnlyPositive = true;
				if (param.mLongSwitch == "only_negative")
					mSemiSupervisedInduceOnlyNegative = true;
				if (param.mShortSwitch == "V")
					mVerbose = true;
				if (param.mShortSwitch == "W")
					mMinimalOutput = true;
				if (param.mLongSwitch == "sequence_token")
					mSequenceToken = true;
				if (param.mLongSwitch == "sequence_multi_line")
					mSequenceMultiLine = true;
				if (param.mLongSwitch == "sequence_pairwise_interaction")
					mSequencePairwiseInteraction = true;
				if (param.mShortSwitch == "X")
					mSharedNeighborhood = true;
				if (param.mLongSwitch == "binary_file_type")
					mBinaryFormat = true;
				if (param.mLongSwitch == "force_approximate")
					mForceApproximate = true;
				if (param.mLongSwitch == "no_neighborhood_cache")
					mNoNeighborhoodCache = true;
				if (param.mLongSwitch == "no_minhash_cache")
					mNoMinHashCache = true;
			}

			if (param.mShortSwitch == "a")
				mAction = param.mValue;
			if (param.mShortSwitch == "A")
				mRowIndexFileName = param.mValue;
			if (param.mShortSwitch == "B")
				mColIndexFileName = param.mValue;
			if (param.mShortSwitch == "m")
				mModelFileName = param.mValue;
			if (param.mShortSwitch == "f")
				mFileType = param.mValue;
			if (param.mShortSwitch == "k")
				mKernelType = param.mValue;
			if (param.mShortSwitch == "g")
				mGraphType = param.mValue;
			if (param.mShortSwitch == "o")
				mOpenBabelFormat = param.mValue;
			if (param.mShortSwitch == "s")
				mSuffix = param.mValue;
			if (param.mShortSwitch == "r")
				mRadius = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "d")
				mDistance = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "v")
				mVertexDegreeThreshold = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "M")
				mSequenceDegree = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "b")
				mHashBitSize = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "l")
				mLambda = stream_cast<double>(param.mValue);
			if (param.mShortSwitch == "e")
				mEpochs = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "c")
				mCrossValidationNumFolds = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "p")
				mNumPoints = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "R")
				mRandomSeed = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "S")
				mSemiSupervisedNumIterations = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "T")
				mSemiSupervisedThreshold = stream_cast<double>(param.mValue);
			if (param.mShortSwitch == "N")
				mLMDSNumRandomRestarts = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "I")
				mLMDSNumIterations = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "D")
				mLMDSDimensionality = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "E")
				mLMDSIterationEpsilon = stream_cast<double>(param.mValue);
			if (param.mShortSwitch == "K")
				mLMDSNeighborhoodSize = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "n")
				mLMDSNonNeighborhoodSize = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "G")
				mLMDSNeighborhoodSizeRange = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "U")
				mLMDSTau = stream_cast<double>(param.mValue);
			if (param.mShortSwitch == "u")
				mLMDSTauExponentRange = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "O")
				mSparsificationNumIterations = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "C")
				mTopologicalRegularizationNumNeighbors = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "L")
				mTopologicalRegularizationRate = stream_cast<double>(param.mValue);
			if (param.mShortSwitch == "F")
				mNumHashFunctions = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "z")
				mMaxSizeBin = stream_cast<double>(param.mValue);
			if (param.mLongSwitch == "eccess_neighbour_size_factor")
				mEccessNeighbourSizeFactor = stream_cast<double>(param.mValue);
			if (param.mShortSwitch == "Q")
				mSampleSize = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "x")
				mNumNearestNeighbors = stream_cast<unsigned>(param.mValue);
			if (param.mLongSwitch == "fraction_center_scan")
				mFractionCenterScan = stream_cast<double>(param.mValue);
			if (param.mShortSwitch == "w")
				mMaxIntersectionSize = stream_cast<unsigned>(param.mValue);
			if (param.mShortSwitch == "Y")
				mNumLineSearchIterations = stream_cast<unsigned>(param.mValue);
			if (param.mLongSwitch == "num_repeat_hash_functions")
				mNumRepeatsHashFunction = stream_cast<unsigned>(param.mValue);
			if (param.mLongSwitch == "cluster_threshold")
				mClusterThreshold = stream_cast<double>(param.mValue);
			if (param.mLongSwitch == "cluster_type")
				mClusterType = param.mValue;
			if (param.mLongSwitch == "semi_supervised_alpha")
				mSemiSupervisedAlpha = stream_cast<double>(param.mValue);
			//DD kernel family additional parameters
			if (param.mLongSwitch == "tree_lambda")
				mTreeLambda = stream_cast<double>(param.mValue);
			if (param.mLongSwitch == "radius_two")
				mRadiusTwo = stream_cast<double>(param.mValue);
		}

		//convert action string to action code
		if (mAction == "")
			mActionCode = NULL_ACTION;
		else if (mAction == "TRAIN")
			mActionCode = TRAIN;
		else if (mAction == "TEST")
			mActionCode = TEST;
		else if (mAction == "CROSS_VALIDATION")
			mActionCode = CROSS_VALIDATION;
		else if (mAction == "BIAS_VARIANCE_DECOMPOSITION")
			mActionCode = BIAS_VARIANCE_DECOMPOSITION;
		else if (mAction == "PARAMETERS_OPTIMIZATION")
			mActionCode = PARAMETERS_OPTIMIZATION;
		else if (mAction == "LEARNING_CURVE")
			mActionCode = LEARNING_CURVE;
		else if (mAction == "TEST_PART")
			mActionCode = TEST_PART;
		else if (mAction == "FEATURE")
			mActionCode = FEATURE;
		else if (mAction == "FEATURE_PART")
			mActionCode = FEATURE_PART;
		else if (mAction == "FEATURE_SCALED")
			mActionCode = FEATURE_SCALED;
		else if (mAction == "MATRIX")
			mActionCode = MATRIX;
		else if (mAction == "EMBED")
			mActionCode = EMBED;
		else if (mAction == "TARGET_ALIGNMENT")
			mActionCode = TARGET_ALIGNMENT;
		else if (mAction == "CLUSTER")
			mActionCode = CLUSTER;
		else if (mAction == "NEAREST_NEIGHBOR")
			mActionCode = NEAREST_NEIGHBOR;
		else if (mAction == "MIN_HASH")
			mActionCode = MIN_HASH;
		else if (mAction == "SEMI_SUPERVISED")
			mActionCode = SEMI_SUPERVISED;
		else
			throw range_error("ERROR2.46: Unrecognized action: <" + mAction + ">");

		//convert file type string to file type code
		if (mFileType == "GRAPH")
			mFileTypeCode = GRAPH;
		else if (mFileType == "SPARSE_VECTOR")
			mFileTypeCode = SPARSE_VECTOR;
#ifdef USEOBABEL
		else if (mFileType == "MOLECULAR_GRAPH")
			mFileTypeCode = MOLECULAR_GRAPH;
#endif
		else if (mFileType == "SEQUENCE")
			mFileTypeCode = SEQUENCE;
		else
			throw range_error("ERROR2.46: Unrecognized file type: <" + mFileType + ">");

		//check for help request
		for (unsigned i = 0; i < options.size(); ++i) {
			if (options[i] == "-h" || options[i] == "--help") {
				Usage(argv[0], "COMPACT");
				exit(1);
			}
			if (options[i] == "-H" || options[i] == "--Help") {
				Usage(argv[0], "EXTENDED");
				exit(1);
			}
		}

		//check that set parameters are compatible
		if (mInputDataFileName == "")
			throw range_error("ERROR2.4: -i <input data file name> is missing.");

		if ((mAction == "TRAIN" || mAction == "CROSS_VALIDATION" || mAction == "LEARNING_CURVE" || mAction == "TARGET_ALIGNMENT") && mTargetFileName == "")
			throw range_error("ERROR2.5: -t <target file name> is missing.");

#ifndef USEOBABEL
		if (mFileType == "MOLECULAR_GRAPH")
		throw range_error("ERROR5: -f MOLECULAR_GRAPH is not enabled if compilation flag USEOBABEL is not set. Please recompile using -DUSEOBABEL");
#endif
	}
};

//------------------------------------------------------------------------------------------------------------------------
class Kernel {
public:
	NSPDK_FeatureGenerator* mpFeatureGenerator;
	Parameters* mpParameters;
public:
	void Init(Parameters* apParameters) {
		mpParameters = apParameters;
		if (mpParameters->mKernelType == "NSPDK")
			mpFeatureGenerator = new ANSPDK_FeatureGenerator("anspdk");
		else if (mpParameters->mKernelType == "PBK")
			mpFeatureGenerator = new PBK_FeatureGenerator("pbk");
		else if (mpParameters->mKernelType == "WDK")
			mpFeatureGenerator = new WDK_FeatureGenerator("wdk");
		else if (mpParameters->mKernelType == "USPK")
			mpFeatureGenerator = new USPK_FeatureGenerator("uspk");
		else if (mpParameters->mKernelType == "DDK")
					mpFeatureGenerator = new DDkernel_FeatureGeneratorNew("ddk");
		else if (mpParameters->mKernelType == "NSDDK")
					mpFeatureGenerator = new NSDDkernel_FeatureGenerator("nsddk");
		else if (mpParameters->mKernelType == "ANSDDK")
					mpFeatureGenerator = new ANSDDkernel_FeatureGenerator("ansddk");
		else if (mpParameters->mKernelType == "SK")
			mpFeatureGenerator = new SK_FeatureGenerator("sk");
		else
			throw range_error("ERROR2.1: Unknown kernel type: " + mpParameters->mKernelType);

		ParametersSetup();
	}

	void ParametersSetup() {
#ifdef DEBUGON
		mpFeatureGenerator->set_flag("verbosity", stream_cast<string>(1));
#endif
		if (mpParameters->mVerbose)
			mpFeatureGenerator->set_flag("verbosity", stream_cast<string>(1));
		if (mpParameters->mMinKernel)
			mpFeatureGenerator->set_flag("min_kernel", "true");
		if (mpParameters->mKernelNoNormalization)
			mpFeatureGenerator->set_flag("normalization", "false");
		mpFeatureGenerator->set_flag("radius", stream_cast<string>(mpParameters->mRadius));
		mpFeatureGenerator->set_flag("distance", stream_cast<string>(mpParameters->mDistance));
		mpFeatureGenerator->set_flag("hash_bit_size", stream_cast<string>(mpParameters->mHashBitSize));
		unsigned bitmask = (2 << (mpParameters->mHashBitSize - 1)) - 1;
		mpFeatureGenerator->set_flag("hash_bit_mask", stream_cast<string>(bitmask));

		//if type of kernel PBK then also perform the following initializations
		if (mpParameters->mKernelType == "PBK") {
			mpFeatureGenerator->set_flag("lower_vertex_degree_threshold", stream_cast<string>(mpParameters->mVertexDegreeThreshold));
			mpFeatureGenerator->set_flag("vertex_degree_threshold", stream_cast<string>(mpParameters->mVertexDegreeThreshold));

			if (mpParameters->mMinKernel)
				dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mWDK.set_flag("min_kernel", "true");
			if (mpParameters->mKernelNoNormalization)
				dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mWDK.set_flag("normalization", "false");
			dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mWDK.set_flag("radius", stream_cast<string>(mpParameters->mRadius));
			dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mWDK.set_flag("distance", stream_cast<string>(mpParameters->mDistance));
			dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mWDK.set_flag("hash_bit_size", stream_cast<string>(mpParameters->mHashBitSize));
			bitmask = (2 << (mpParameters->mHashBitSize - 1)) - 1;
			dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mWDK.set_flag("hash_bit_mask", stream_cast<string>(bitmask));
			dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mWDK.set_flag("lower_vertex_degree_threshold", stream_cast<string>(mpParameters->mVertexDegreeThreshold));

			if (mpParameters->mMinKernel)
				dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mANSPDK.set_flag("min_kernel", "true");
			if (mpParameters->mKernelNoNormalization)
				dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mANSPDK.set_flag("normalization", "false");
			dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mANSPDK.set_flag("radius", stream_cast<string>(mpParameters->mRadius));
			dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mANSPDK.set_flag("distance", stream_cast<string>(mpParameters->mDistance));
			dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mANSPDK.set_flag("hash_bit_size", stream_cast<string>(mpParameters->mHashBitSize));
			bitmask = (2 << (mpParameters->mHashBitSize - 1)) - 1;
			dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mANSPDK.set_flag("hash_bit_mask", stream_cast<string>(bitmask));
			dynamic_cast<PBK_FeatureGenerator*>(mpFeatureGenerator)->mANSPDK.set_flag("vertex_degree_threshold", stream_cast<string>(mpParameters->mVertexDegreeThreshold));

		}

		//DD kernel family additional parameters
		if (mpParameters->mKernelType == "DDK"){
			dynamic_cast<DDkernel_FeatureGeneratorNew*>(mpFeatureGenerator)->set_flag("mTreeLambda", stream_cast<string>(mpParameters->mTreeLambda));
		}
		if (mpParameters->mKernelType == "NSDDK" || mpParameters->mKernelType == "ANSDDK"){
			dynamic_cast<NSDDkernel_FeatureGenerator*>(mpFeatureGenerator)->set_flag("mTreeLambda", stream_cast<string>(mpParameters->mTreeLambda));
			dynamic_cast<NSDDkernel_FeatureGenerator*>(mpFeatureGenerator)->set_flag("mRadiusTwo", stream_cast<string>(mpParameters->mRadiusTwo));

		}

		/*		if (!mpParameters->mMinimalOutput) {
		 cout << SEP << endl << "Kernel parameters" << endl << SEP << endl;
		 mpFeatureGenerator->OutputParameters(cout);
		 cout << SEP << endl;
		 }
		 */
	}

	void GenerateFeatureVector(GraphClass& aG, SVector& oX) {
		mpFeatureGenerator->generate_feature_vector(aG, oX);
	}

	void GenerateVertexFeatureVector(GraphClass& aG, vector<SVector>& oXList) {
		mpFeatureGenerator->generate_vertex_feature_vector(aG, oXList);
	}

	double ComputeKernel(GraphClass& aG, GraphClass& aM) {
		SVector xg;
		SVector xm;
		GenerateFeatureVector(aG, xg);
		GenerateFeatureVector(aM, xm);
		return ComputeKernel(xg, xm);
	}

	double ComputeKernel(const SVector& aX, const SVector& aZ) {
		return dot(aX, aZ);
	}
};

//------------------------------------------------------------------------------------------------------------------------
class Data {
public:
	Parameters* mpParameters;
	Kernel mKernel;
	vector<double> mTargetList;
	vector<SVector> mVectorList;
	map<unsigned, SVector> mFeatureCorrelationMatrix;
	vector<unsigned> mRowIndexList;
	vector<unsigned> mColIndexList;
public:
	Data() {
	}

	void Init(Parameters* apParameters) {
		mpParameters = apParameters;
		mKernel.Init(mpParameters);
	}

	double ComputeKernel(unsigned i, unsigned j) {
		if (i > Size() || j > Size())
			throw range_error("ERROR3.1: Kernel computation for instances out of range. Available range: 0-" + stream_cast<string>(Size()) + " but asked for kernel between instances with index " + stream_cast<string>(i) + "," + stream_cast<string>(j));
		return mKernel.ComputeKernel(mVectorList[i], mVectorList[j]);
	}

	void ComputeFeatureCorrelationMatrix() {
		mFeatureCorrelationMatrix.clear();
		igzstream fin;
		fin.open(mpParameters->mInputDataFileName.c_str());
		if (!fin)
			throw range_error("ERROR2.11: Cannot open file: " + mpParameters->mInputDataFileName);

		{
			if (!mpParameters->mMinimalOutput)
				cout << endl << "Processing file:" << mpParameters->mInputDataFileName << " for feature correlation" << endl;
			ProgressBar pb;
			bool valid_input = true;
			while (!fin.eof() && valid_input) {
				GraphClass g;
				SetGraphFromFile(fin, g);
				if (!g.IsEmpty()) {
					if (!mpParameters->mMinimalOutput)
						pb.Count();
					//extract from a graph the feature vector for each vertex (manage directed graph)
					vector<SVector> graph_vertex_vector_list;
					mKernel.GenerateVertexFeatureVector(g, graph_vertex_vector_list);
					unsigned size = mpParameters->mGraphType == "DIRECTED" ? graph_vertex_vector_list.size() / 2 : graph_vertex_vector_list.size();
					for (unsigned vertex_id = 0; vertex_id < size; ++vertex_id) {
						SVector x = graph_vertex_vector_list[vertex_id];
						if (mpParameters->mGraphType == "DIRECTED")
							x.add(graph_vertex_vector_list[vertex_id + size]);
					}

					//impute to each feature the occurrences of the other features associated to the same single vertex
					//for each vertex
					for (unsigned i = 0; i < graph_vertex_vector_list.size(); ++i) {
						SVector& x = graph_vertex_vector_list[i];
						vector<pair<int, double> > vec = x.unpack();
						//for each feature
						for (unsigned j = 0; j < vec.size(); ++j) {
							unsigned hash_code = vec[j].first;
							if (mFeatureCorrelationMatrix.count(hash_code) == 0)
								mFeatureCorrelationMatrix.insert(make_pair(hash_code, x));
							else
								mFeatureCorrelationMatrix[hash_code].add(x);
						}
					}
				} else
					valid_input = false;
			}
		}

		//row normalize correlation matrix
		for (map<unsigned, SVector>::iterator it = mFeatureCorrelationMatrix.begin(); it != mFeatureCorrelationMatrix.end(); ++it) {
			unsigned id = it->first;
			SVector& x = it->second;

			//extract threshold for k largest correlation
			vector<pair<int, double> > vec = x.unpack();
			vector<double> corr;
			for (unsigned j = 0; j < vec.size(); ++j)
				corr.push_back(-vec[j].second);
			sort(corr.begin(), corr.end());
			double threshold = -corr[mpParameters->mTopologicalRegularizationNumNeighbors];
			//trim vector to the top mpParameters->mTopologicalRegularizationNumNeighbors
			SVector x_new;
			for (unsigned j = 0; j < vec.size(); ++j)
				if (vec[j].second > threshold)
					x_new.set(vec[j].first, vec[j].second);
			//retain only the top mpParameters->mTopologicalRegularizationNumNeighbors related features by replacing the vector row
			x = x_new;

			//create the normalized Laplacian
			x.set(id, 0); //remove counts to self
			x.normalize();
			x.scale(-1);
			x.set(id, 1); //set counts to self to 1
		}

	}

	void LoadIndex() {
		if (mpParameters->mRowIndexFileName != "")
			if (mRowIndexList.size() == 0)
				LoadUnsignedList(mpParameters->mRowIndexFileName, mRowIndexList);
		if (mpParameters->mColIndexFileName != "")
			if (mColIndexList.size() == 0)
				LoadUnsignedList(mpParameters->mColIndexFileName, mColIndexList);
	}

	void LoadTarget() {
		if (mTargetList.size() == 0)
			LoadRealList(mpParameters->mTargetFileName, mTargetList);
	}

	void LoadRealList(string aFileName, vector<double>& oList) {
		oList.clear();
		if (!mpParameters->mMinimalOutput)
			cout << endl << "Reading file: " << aFileName << " ..";
		ifstream fin;
		fin.open(aFileName.c_str());
		if (!fin)
			throw range_error("ERROR2.235: Cannot open file:" + aFileName);
		while (!fin.eof()) {
			string line;
			getline(fin, line);
			stringstream ss;
			ss << line << endl;
			while (!ss.eof()) {
				double value(0);
				ss >> value;
				if (ss.good()) {
					oList.push_back(value);
				}
			}
		}
		fin.close();
		if (!mpParameters->mMinimalOutput)
			cout << ".. read: " << oList.size() << " values." << endl;
	}

	void LoadUnsignedList(string aFileName, vector<unsigned>& oList) {
		oList.clear();
		if (!mpParameters->mMinimalOutput)
			cout << endl << "Reading file: " << aFileName << " ..";
		ifstream fin;
		fin.open(aFileName.c_str());
		if (!fin)
			throw range_error("ERROR2.236: Cannot open file:" + aFileName);
		while (!fin.eof()) {
			string line;
			getline(fin, line);
			stringstream ss;
			ss << line << endl;
			while (!ss.eof()) {
				unsigned value(0);
				ss >> value;
				if (ss.good()) {
					oList.push_back(value);
				}
			}
		}
		fin.close();
		if (!mpParameters->mMinimalOutput)
			cout << ".. read: " << oList.size() << " values." << endl;
	}

	void LoadGspanList(vector<GraphClass>& aGraphList) {
		mVectorList.clear();
		for (unsigned i = 0; i < aGraphList.size(); ++i) {
			SVector x;
			mKernel.GenerateFeatureVector(aGraphList[i], x);
			mVectorList.push_back(x);
		}
	}

	void LoadData() {
		mVectorList.clear();
		mKernel.ParametersSetup();
		igzstream fin;
		fin.open(mpParameters->mInputDataFileName.c_str());
		if (!fin)
			throw range_error("ERROR2.11: Cannot open file: " + mpParameters->mInputDataFileName);
		ProgressBar pb;
		if (!mpParameters->mMinimalOutput)
			cout << "Processing file: " << mpParameters->mInputDataFileName << endl;
		bool valid_input = true;
		while (!fin.eof() && valid_input) {
			switch (mpParameters->mFileTypeCode) {
			case GRAPH:
#ifdef USEOBABEL
			case MOLECULAR_GRAPH:
#endif
				case SEQUENCE: {
				vector<GraphClass> g_list(BUFFER_SIZE);
				unsigned i = 0;
				while (i < BUFFER_SIZE && !fin.eof() && valid_input) {
					SetGraphFromFile(fin, g_list[i]);
					if (g_list[i].IsEmpty())
						valid_input = false;
					else
						i++;
				}
#pragma omp parallel for schedule(dynamic,1) ordered
				for (unsigned j = 0; j < i; j++) {
#pragma omp ordered
					{
						SVector x;
						mKernel.GenerateFeatureVector(g_list[j], x);
						mVectorList.push_back(x);
						pb.Count();
					}
				}
			}
				break;
			case SPARSE_VECTOR: {
				SVector x;
				if (mpParameters->mBinaryFormat)
					SetVectorFromSparseVectorBinaryFile(fin, x);
				else
					SetVectorFromSparseVectorAsciiFile(fin, x);
				if (x.size() > 0) {
					mVectorList.push_back(x);
					if (!mpParameters->mMinimalOutput)
						pb.Count();
				}
			}
				break;
			default:
				throw range_error("ERROR2.45: file type not recognized: " + mpParameters->mFileType);
			}
		}

		//additional files to load if required
		if (mpParameters->mRowIndexFileName != "" || mpParameters->mColIndexFileName != "")
			LoadIndex();

		if (mpParameters->mTargetFileName != "")
			LoadTarget();

		//additional operations to perform on data if required
		if (mpParameters->mTopologicalRegularizationNumNeighbors != 0)
			ComputeFeatureCorrelationMatrix();

		if (mRowIndexList.size() == 0) {
			if (!mpParameters->mMinimalOutput)
				cout << endl << "No row index list specified. Assuming all " << Size() << " row indices as valid." << endl;
			for (unsigned i = 0; i < Size(); ++i)
				mRowIndexList.push_back(i);
		}
		if (mColIndexList.size() == 0) {
			if (!mpParameters->mMinimalOutput)
				cout << endl << "No col index list specified. Assuming all " << Size() << " col indices as valid." << endl;
			for (unsigned i = 0; i < Size(); ++i)
				mColIndexList.push_back(i);
		}
	}

	void SetGraphFromFile(istream& in, GraphClass& oG) {
		switch (mpParameters->mFileTypeCode) {
		case GRAPH:
			SetGraphFromGraphGspanFile(in, oG);
			break;
#ifdef USEOBABEL
		case MOLECULAR_GRAPH:
			SetGraphFromGraphOpenBabelFile(in, oG, mpParameters->mOpenBabelFormat);
			break;
#endif
		case SEQUENCE: {
			if (mpParameters->mSequenceToken && mpParameters->mSequenceMultiLine)
				SetGraphFromSequenceMultiLineTokenFile(in, oG);
			if (mpParameters->mSequenceToken && !mpParameters->mSequenceMultiLine)
				SetGraphFromSequenceTokenFile(in, oG);
			if (!mpParameters->mSequenceToken && mpParameters->mSequenceMultiLine)
				SetGraphFromSequenceMultiLineFile(in, oG);
			if (!mpParameters->mSequenceToken && !mpParameters->mSequenceMultiLine)
				SetGraphFromSequenceFile(in, oG);
			break;
		}
		default:
			throw range_error("ERROR2.45: file type not recognized: " + mpParameters->mFileType);
		}
		//HACK**************************************************************************************************
		/*static int counter = -1;
		 if (oG.VertexSize() > 0) {
		 counter++;
		 string gid = stream_cast<string>(counter);
		 vector<unsigned> view_point_list;
		 oG.ComputePairwiseDistanceInformation(mpParameters->mRadius + mpParameters->mDistance);
		 oG.SaveAsMatlabFile(gid, view_point_list);
		 }
		 */
		//******************************************************************************************************
	}

	void SetGraphFromGraphOpenBabelFile(istream& in, GraphClass& oG, string aFormat) {
		static OpenBabelConverter molecule_converter;
		if (in.good() && !in.eof()) {
			//molecule_converter.ConvertOpenBabelFormatToGraph(&in, oG, aFormat); //OLD version
			molecule_converter.ConvertOpenBabelFormatToGraph(mpParameters->mInputDataFileName, oG, aFormat);
		}
	}

	void SetVectorFromSparseVectorAsciiFile(istream& in, SVector& aX) {
		string line;
		getline(in, line);
		if (line == "")
			return;
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

	void SetVectorFromSparseVectorBinaryFile(istream& in, SVector& aX) {
		aX.load(in);
	}

	void SetGraphFromSequenceFile(istream& in, GraphClass& oG) {
		vector<vector<unsigned> > vertex_component_list;
		vector<unsigned> vertex_component;
		bool graph_disconnect = true;
		string line;
		getline(in, line);
		if (line == "")
			return;
		unsigned vertex_counter = 0;
		for (unsigned position_counter = 0; position_counter < line.length(); position_counter++) {
			char char_label = line.at(position_counter);
			string label = stream_cast<string>(char_label);
			if (label == "|") {
				graph_disconnect = true;
				vertex_component_list.push_back(vertex_component);
				vertex_component.clear();
			} else {
				vertex_component.push_back(vertex_counter);
				AddVertexAndEdgesForSequence(oG, label, vertex_counter, graph_disconnect);
				vertex_counter++;
				graph_disconnect = false;
			}
		}
		vertex_component_list.push_back(vertex_component);
		ManageDirectedAndMultiComponent(oG, vertex_component_list);
	}

	void SetGraphFromSequenceTokenFile(istream& in, GraphClass& oG) {
		vector<vector<unsigned> > vertex_component_list;
		vector<unsigned> vertex_component;
		bool graph_disconnect = true;
		string line;
		getline(in, line);
		if (line == "")
			return;
		stringstream ss;
		ss << line << endl;
		unsigned vertex_counter = 0;
		while (!ss.eof() && ss.good()) {
			//add vertex
			string label;
			ss >> label;
			if (label != "") {
				if (label == "|") {
					graph_disconnect = true;
					vertex_component_list.push_back(vertex_component);
					vertex_component.clear();
				} else {
					vertex_component.push_back(vertex_counter);
					AddVertexAndEdgesForSequence(oG, label, vertex_counter, graph_disconnect);
					vertex_counter++;
					graph_disconnect = false;
				}
			}
		}
		ManageDirectedAndMultiComponent(oG, vertex_component_list);
	}

	void SetGraphFromSequenceMultiLineFile(istream& in, GraphClass& oG) {
		vector<vector<unsigned> > vertex_component_list;
		vector<unsigned> vertex_component;
		bool graph_disconnect = true;
		vector<string> line(mpParameters->mSequenceDegree);
		for (unsigned i = 0; i < mpParameters->mSequenceDegree; ++i) {
			getline(in, line[i]);
		}
		unsigned sequence_length = line[0].length();
		unsigned vertex_counter = 0;
		for (unsigned j = 0; j < sequence_length; j++) {
			for (unsigned i = 0; i < mpParameters->mSequenceDegree; ++i) {
				if (line[i] == "")
					return;
				char char_label = line[i].at(j);
				string label = stream_cast<string>(char_label);
				if (label == "|") {
					graph_disconnect = true;
					vertex_component_list.push_back(vertex_component);
					vertex_component.clear();
				} else {
					vertex_component.push_back(vertex_counter);
					AddVertexAndEdgesForSequence(oG, label, vertex_counter, graph_disconnect);
					vertex_counter++;
					graph_disconnect = false;
				}
			}
		}
		ManageDirectedAndMultiComponent(oG, vertex_component_list);
	}

	void SetGraphFromSequenceMultiLineTokenFile(istream& in, GraphClass& oG) {
		vector<vector<unsigned> > vertex_component_list;
		vector<unsigned> vertex_component;
		bool graph_disconnect = true;
		vector<vector<string> > multi_line(mpParameters->mSequenceDegree);
		for (unsigned i = 0; i < mpParameters->mSequenceDegree; ++i) {
			vector<string> single_line;
			string line;
			getline(in, line);
			stringstream ss;
			ss << line << endl;
			while (!ss.eof() && ss.good()) {
				string label;
				ss >> label;
				if (label != "")
					single_line.push_back(label);
			}
			multi_line.push_back(single_line);
		}

		unsigned sequence_length = multi_line[0].size();
		unsigned vertex_counter = 0;
		for (unsigned j = 0; j < sequence_length; j++) {
			for (unsigned i = 0; i < mpParameters->mSequenceDegree; ++i) {
				string label = multi_line[i][j];
				if (label == "|") {
					graph_disconnect = true;
					vertex_component_list.push_back(vertex_component);
					vertex_component.clear();
				} else {
					vertex_component_list.push_back(vertex_component);
					AddVertexAndEdgesForSequence(oG, label, vertex_counter, graph_disconnect);
					vertex_counter++;
					graph_disconnect = false;
				}
			}
		}
		ManageDirectedAndMultiComponent(oG, vertex_component_list);
	}

	void ManageDirectedAndMultiComponent(GraphClass& oG, vector<vector<unsigned> >& aVertexComponentList) {
		if (mpParameters->mGraphType == "DIRECTED" && aVertexComponentList.size() > 1) {
			//create the components for the reversed direction graph
			vector<vector<unsigned> > reverse_vertex_component_list;
			for (unsigned u = 0; u < aVertexComponentList.size(); ++u) {
				vector<unsigned> reverse_vertex_component;
				for (unsigned t = 0; t < aVertexComponentList[u].size(); ++t) {
					unsigned id = aVertexComponentList[u][t];
					reverse_vertex_component.push_back(id + oG.VertexSize());
				}
				reverse_vertex_component_list.push_back(reverse_vertex_component);
			}
			//add the reverse direction components to the component list
			aVertexComponentList.insert(aVertexComponentList.end(), reverse_vertex_component_list.begin(), reverse_vertex_component_list.end());
		}

		if (mpParameters->mGraphType == "DIRECTED")
			AddReverseGraph(oG);

		if (mpParameters->mSequencePairwiseInteraction)
			if (aVertexComponentList.size() > 1)
				AddAbstractConnections(oG, aVertexComponentList);
	}

	void AddVertexAndEdgesForSequence(GraphClass& oG, string aLabel, unsigned aVertexCounter, bool aGraphDisconnect) {
		//set (once) boolean status vectors
		static vector<bool> vertex_status(5, false);
		vertex_status[0] = true; //kernel point
		vertex_status[1] = true; //kind
		vertex_status[2] = true; //viewpoint
		//vertex_status[3] = false; //dead
		//vertex_status[4] = false; //abstraction

		static vector<bool> edge_status(3, false);
		//edge_status[0] = false; //edge dead
		//edge_status[1] = false; //edge abstraction_of
		//edge_status[2] = false; //edge part_of

		unsigned real_vertex_index = oG.InsertVertex();
		vector<string> vertex_symbolic_attribute_list;
		vertex_symbolic_attribute_list.push_back(aLabel);
		oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
		oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);
		if (aVertexCounter % mpParameters->mSequenceDegree != 0) {
			oG.SetVertexViewPoint(real_vertex_index, false);
		}
		if (aGraphDisconnect == false) {
			//add edge
			unsigned real_src_index;
			unsigned real_dest_index;
			if (aVertexCounter % mpParameters->mSequenceDegree != 0) { //add edge to preceding vertex
				real_src_index = aVertexCounter;
				real_dest_index = aVertexCounter - 1;
			} else { //vertices with position multiple than mSequenceDegree are connected in sequence
				real_src_index = aVertexCounter;
				real_dest_index = aVertexCounter - mpParameters->mSequenceDegree;
			}
			vector<string> edge_symbolic_attribute_list;
			edge_symbolic_attribute_list.push_back("-"); //NOTE: default edge label is '-'
			unsigned edge_index = oG.InsertEdge(real_src_index, real_dest_index);
			oG.SetEdgeSymbolicAttributeList(edge_index, edge_symbolic_attribute_list);
			oG.SetEdgeStatusAttributeList(edge_index, edge_status);
			if (mpParameters->mGraphType == "UNDIRECTED" || aVertexCounter % mpParameters->mSequenceDegree != 0) { //add reverse edge in case the graph is undirected or for the attribute vertices, so that in both directions these are accessible
				unsigned reverse_edge_index = oG.InsertEdge(real_dest_index, real_src_index);
				oG.SetEdgeSymbolicAttributeList(reverse_edge_index, oG.GetEdgeSymbolicAttributeList(edge_index));
				oG.SetEdgeStatusAttributeList(reverse_edge_index, oG.GetEdgeStatusAttributeList(edge_index));
			}
		}
	}

	void AddAbstractConnections(GraphClass& oG, vector<vector<unsigned> >& aVertexComponentList) {
		//set (once) boolean status vectors
		static vector<bool> vertex_status(5, false);
		//vertex_status[0] = false; //kernel point
		//vertex_status[1] = false; //kind
		//vertex_status[2] = false; //viewpoint
		//vertex_status[3] = false; //dead
		vertex_status[4] = true; //abstraction

		static vector<bool> edge_status(3, false);
		//edge_status[0] = false; //edge dead
		//edge_status[1] = false; //edge abstraction_of
		//edge_status[2] = false; //edge part_of

		//for all pairs of components
		for (unsigned i = 0; i < aVertexComponentList.size(); ++i) {
			for (unsigned j = i + 1; j < aVertexComponentList.size(); ++j) {
				//join all vertices in one component to all other vertices in the other component
				//add 1 abstract vertex
				unsigned real_vertex_index = oG.InsertVertex();
				vector<string> vertex_symbolic_attribute_list;
				vertex_symbolic_attribute_list.push_back("^L");
				oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
				oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);

				for (unsigned ii = 0; ii < aVertexComponentList[i].size(); ii++) { //add part_of edges
					unsigned real_src_index = real_vertex_index;
					unsigned real_dest_index = aVertexComponentList[i][ii];
					//only add edges to point vertices
					if (oG.GetVertexViewPoint(real_dest_index)) {
						vector<string> edge_symbolic_attribute_list;
						edge_symbolic_attribute_list.push_back("@-");
						unsigned edge_index = oG.InsertEdge(real_src_index, real_dest_index);
						oG.SetEdgeSymbolicAttributeList(edge_index, edge_symbolic_attribute_list);
						oG.SetEdgeStatusAttributeList(edge_index, edge_status);
						oG.SetEdgePartOf(edge_index, true);
					}
				}

				for (unsigned jj = 0; jj < aVertexComponentList[j].size(); jj++) { //add abstraction_of edges
					unsigned real_src_index = real_vertex_index;
					unsigned real_dest_index = aVertexComponentList[j][jj];
					if (oG.GetVertexViewPoint(real_dest_index)) {
						vector<string> edge_symbolic_attribute_list;
						edge_symbolic_attribute_list.push_back("^-");
						unsigned edge_index = oG.InsertEdge(real_src_index, real_dest_index);
						oG.SetEdgeSymbolicAttributeList(edge_index, edge_symbolic_attribute_list);
						oG.SetEdgeStatusAttributeList(edge_index, edge_status);
						oG.SetEdgeAbstractionOf(edge_index, true);
					}
				}
			}
		}
	}

	void SetGraphFromGraphGspanFile(istream& in, GraphClass& oG) {
		//status
		vector<bool> vertex_status(5, false);
		vertex_status[0] = true; //kernel point
		vertex_status[1] = true; //kind
		vertex_status[2] = true; //viewpoint
		vertex_status[3] = false; //dead
		vertex_status[4] = false; //abstraction

		vector<bool> edge_status(3, false);
		edge_status[0] = false; //edge dead
		edge_status[1] = false; //edge abstraction_of
		edge_status[2] = false; //edge part_of

		map<string, int> index_map_nominal_to_real;
		string gid;
		string line;
		unsigned line_counter = 0;
		bool instance_started = false;
		do {
			char c = in.peek();
			if (c == 't' && instance_started == true)
				break;

			line_counter++;
			getline(in, line);
			if (line == "")
				break;
			stringstream ss;
			ss << line << endl;
			char code;
			ss >> code;
			if (code == 't') {
				instance_started = true;
			} else if (code == 'v' || code == 'V' || code == 'W') {
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
				if (code == 'V')
					oG.SetVertexViewPoint(real_vertex_index, false);
				if (code == 'W') {
					oG.SetVertexViewPoint(real_vertex_index, false);
					oG.SetVertexKernelPoint(real_vertex_index, false);
				}
				char vertex_abstraction_code = label.at(0);
				if (vertex_abstraction_code == '^') {
					oG.SetVertexAbstraction(real_vertex_index, true);
					oG.SetVertexViewPoint(real_vertex_index, false);
					oG.SetVertexKernelPoint(real_vertex_index, false);
				}
			} else if (code == 'e') {
				//extract src and dest vertex id
				string nominal_src_index, nominal_dest_index;
				string label;
				ss >> nominal_src_index >> nominal_dest_index >> label;
				assert(index_map_nominal_to_real.count(nominal_src_index) > 0);
				assert(index_map_nominal_to_real.count(nominal_dest_index) > 0);
				vector<string> edge_symbolic_attribute_list;
				edge_symbolic_attribute_list.push_back(label);
				unsigned real_src_index = index_map_nominal_to_real[nominal_src_index];
				unsigned real_dest_index = index_map_nominal_to_real[nominal_dest_index];
				unsigned edge_index = oG.InsertEdge(real_src_index, real_dest_index);
				oG.SetEdgeSymbolicAttributeList(edge_index, edge_symbolic_attribute_list);
				oG.SetEdgeStatusAttributeList(edge_index, edge_status);

				char edge_abstraction_code = label.at(0);
				if (edge_abstraction_code == '^')
					oG.SetEdgeAbstractionOf(edge_index, true);
				if (edge_abstraction_code == '@')
					oG.SetEdgePartOf(edge_index, true);

				if (mpParameters->mGraphType == "UNDIRECTED" || edge_abstraction_code == '^' || edge_abstraction_code == '@') { //NOTE: edges that are part of the abstraction mechanism should be treated as undirected
					unsigned reverse_edge_index = oG.InsertEdge(real_dest_index, real_src_index);
					oG.SetEdgeSymbolicAttributeList(reverse_edge_index, oG.GetEdgeSymbolicAttributeList(edge_index));
					oG.SetEdgeStatusAttributeList(reverse_edge_index, oG.GetEdgeStatusAttributeList(edge_index));
				}
			} else {
			} //NOTE: ignore other markers
		} while (!in.eof() && in.good());
		if (mpParameters->mGraphType == "DIRECTED")
			AddReverseGraph(oG);
	}

	void AddReverseGraph(GraphClass& oG) {
		unsigned vsize = oG.VertexSize();
		//add a copy of all vertices
		for (unsigned i = 0; i < vsize; i++) {
			unsigned real_vertex_index = oG.InsertVertex();
			assert(real_vertex_index == i + vsize);
			vector<string> vertex_symbolic_attribute_list = oG.GetVertexSymbolicAttributeList(i);
			for (unsigned t = 0; t < vertex_symbolic_attribute_list.size(); t++) //prepend a prefix to mark the reverse direction
				vertex_symbolic_attribute_list[t] = "r." + vertex_symbolic_attribute_list[t];
			oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
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

	unsigned Size() {
		return mVectorList.size();
	}
};

//------------------------------------------------------------------------------------------------------------------------
/**
 Encapsulates a linear SVM model trainable with stochastic gradient
 descent over graph instances explicitly mapped by the NSPDK kernel
 */
class StochasticGradientDescentSupportVectorMachine {
	/**
	 Data structure to: 1) facilitate rebalancing of dataset by copying
	 multiple times a reference to the instance; and 2) retrieve
	 prediction efficiently by overwriting a reference to the margin
	 list cell element.
	 */
	struct TrainItem {
		int mInstanceID;
		double mTarget;
		double* mpMargin;
		SVector* mpInstance;
	};
public:
	Parameters* mpParameters;
	Data* mpData;
protected:
	double mWScale;
	double mBias;
	SVector mW;

public:
	StochasticGradientDescentSupportVectorMachine() {
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		Clear();
	}

	void Clear() {
		mWScale = 1;
		mBias = 0;
		mW = SVector();
	}

	void VectorElementwiseProductWithModel(SVector& oX) {
		SVector x = elementwise_product(oX, mW);
		x.normalize();
		oX = x;
	}

	/**
	 Prints several informative measures given a list of predictions and a list of true targets
	 */
	void OutputPerformanceMeasures(ostream& out, const vector<double>& aMarginList, const vector<double>& aTargetList) {
		assert(aMarginList.size() == aTargetList.size());
		unsigned size = aMarginList.size();
		unsigned error = 0;
		unsigned correct = 0;
		unsigned tp, tn, fp, fn;
		tp = tn = fp = fn = 0;
		for (unsigned i = 0; i < aMarginList.size(); ++i) {
			double margin = aMarginList[i];
			double prediction = margin > 0 ? 1 : -1;
			double target = aTargetList[i];
			if (prediction != target)
				error++;
			if (prediction == target)
				correct++;
			if (prediction > 0 && target > 0)
				tp++;
			if (prediction > 0 && target < 0)
				fp++;
			if (prediction < 0 && target > 0)
				fn++;
			if (prediction < 0 && target < 0)
				tn++;
		}

		double pprecision = (double) tp / (tp + fp);
		double precall = (double) tp / (tp + fn);
		double pfmeasure = 2 * pprecision * precall / (pprecision + precall);

		double nprecision = (double) tn / (tn + fn);
		double nrecall = (double) tn / (tn + fp);
		double nfmeasure = 2 * nprecision * nrecall / (nprecision + nrecall);

		double bprecision = (pprecision + nprecision) / 2;
		double brecall = (precall + nrecall) / 2;
		double bfmeasure = (pfmeasure + nfmeasure) / 2;

		out << TAB << "Size: " << size << endl;
		out << TAB << "Correct: " << correct << " ( " << correct * 100 / (double) (correct + error) << " %)" << endl;
		out << TAB << "Error: " << error << " ( " << error * 100 / (double) (correct + error) << " %)" << endl;
		out << TAB << "Confusion table:" << endl;
		out << TAB << "TP:" << tp << " FP:" << fp << endl;
		out << TAB << "FN:" << fn << " TN:" << tn << endl;
		out << TAB << "+Precision:" << pprecision << " +Recall:" << precall << " +F-measure:" << pfmeasure << endl;
		out << TAB << "-Precision:" << nprecision << " -Recall:" << nrecall << " -F-measure:" << nfmeasure << endl;
		out << TAB << "bPrecision:" << bprecision << " bRecall:" << brecall << " bF-measure:" << bfmeasure << endl;
	}

	double ComputeBalancedFMeasure(const vector<double>& aMarginList, const vector<double>& aTargetList) {
		assert(aMarginList.size() == aTargetList.size());
		unsigned error = 0;
		unsigned correct = 0;
		unsigned tp, tn, fp, fn;
		tp = tn = fp = fn = 0;
		for (unsigned i = 0; i < aMarginList.size(); ++i) {
			double margin = aMarginList[i];
			double prediction = margin > 0 ? 1 : -1;
			double target = aTargetList[i];
			if (prediction != target)
				error++;
			if (prediction == target)
				correct++;
			if (prediction > 0 && target > 0)
				tp++;
			if (prediction > 0 && target < 0)
				fp++;
			if (prediction < 0 && target > 0)
				fn++;
			if (prediction < 0 && target < 0)
				tn++;
		}

		double pprecision = (double) tp / (tp + fp);
		double precall = (double) tp / (tp + fn);
		double pfmeasure = 2 * pprecision * precall / (pprecision + precall);

		double nprecision = (double) tn / (tn + fn);
		double nrecall = (double) tn / (tn + fp);
		double nfmeasure = 2 * nprecision * nrecall / (nprecision + nrecall);

		double bfmeasure = (pfmeasure + nfmeasure) / 2;

		return bfmeasure;
	}

	void Save(ostream& out) {
		out << "bias " << mBias << endl;
		out << "wscale " << mWScale << endl;
		out << "w " << mW;
	}

	void Save(string aLocalSuffix = string()) {
		string filename = mpParameters->mModelFileName + aLocalSuffix + mpParameters->mSuffix;
		ofstream ofs;
		ofs.open(filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.22: Cannot open file:" + filename);
		if (!mpParameters->mMinimalOutput)
			cout << endl << "Saving model file: " << filename << endl;
		Save(ofs);
	}

	void Load(istream& in) {
		string attribute = "";
		string expected = "";
		in >> attribute >> mBias;
		expected = "bias";
		if (attribute != expected)
			throw range_error("ERROR2.17: Format error: expecting [" + expected + "] but found [" + attribute + "]");
		in >> attribute >> mWScale;
		expected = "wscale";
		if (attribute != expected)
			throw range_error("ERROR2.18: Format error: expecting [" + expected + "] but found [" + attribute + "]");
		in >> attribute >> mW;
		assert(attribute == "w");
	}

	vector<double> Train(vector<double>& aTargetList, vector<unsigned>& aTrainsetIDList) {
		if (aTrainsetIDList.size() != aTargetList.size())
			throw range_error("ERROR2.19: Data list and Target list have not the same size: #data:" + stream_cast<string>(aTrainsetIDList.size()) + " #targets:" + stream_cast<string>(aTargetList.size()));

		vector<SVector*> sv_data_list;
		for (unsigned i = 0; i < aTrainsetIDList.size(); ++i) {
			unsigned id = aTrainsetIDList[i];
			sv_data_list.push_back(&mpData->mVectorList[id]);
		}

		//allocate local target list and compute positive/negative target counts
		vector<double> target_list;
		unsigned p, n;
		p = n = 0;
		for (unsigned i = 0; i < aTargetList.size(); ++i) {
			if (aTargetList[i] > 0)
				p++;
			else
				n++;
			target_list.push_back(aTargetList[i]);
		}

		//if no instance has negative class then generate negative instances with opposite features wrt positive instances
		vector<SVector> synth_neg_sv_data_list(aTrainsetIDList.size());
		if (n == 0) {
			if (!mpParameters->mMinimalOutput)
				cout << "No negative instances: proceeding to generate " << aTrainsetIDList.size() << " negative instances with opposite features wrt positive instances" << endl;
			for (unsigned i = 0; i < aTrainsetIDList.size(); ++i) {
				SVector x = *sv_data_list[i];
				x.scale(-1);
				synth_neg_sv_data_list[i] = x;
				sv_data_list.push_back(&(synth_neg_sv_data_list[i]));
				target_list.push_back(-1);
			}
		}

		//clear margin list and allocate memory
		vector<double> margin_list(target_list.size());
		//...rebalance classes use pointers array to scramble and oversample
		vector<TrainItem> balanced_dataset;
		BalanceDataset(aTrainsetIDList, target_list, margin_list, sv_data_list, balanced_dataset);

		//train on balanced train data
		CoreTrain(balanced_dataset);

		if (!mpParameters->mMinimalOutput) {
			//output statistics on original train data (no class balance)
			OutputModelInfo();
			cout << "Performance on train set:" << endl;
			vector<double> train_margin_list=Test(aTrainsetIDList);
			OutputPerformanceMeasures(cout, train_margin_list, target_list);
		}
		if (n == 0) { //if no instance has negative class then re-test model only on real training data to extract margins and predictions
			margin_list = Test(aTrainsetIDList);
		}
		return margin_list;
	}

	void BalanceDataset(vector<unsigned>& aDatasetIDList, vector<double>& aTargetList, vector<double>& oMarginList, vector<SVector*>& aSVDataList, vector<TrainItem>& oDataset) {
		//compute class distribution
		unsigned p, n;
		p = n = 0;
		for (unsigned i = 0; i < aTargetList.size(); ++i)
			if (aTargetList[i] > 0)
				p++;
			else
				n++;
		if (!mpParameters->mMinimalOutput)
			cout << "Class distribution: " << p + n << " (+:" << p << " -:" << n << ") " << "[+:" << (double) p / (p + n) << " -:" << (double) n / (p + n) << "]" << endl;

		//separate positive from negative instances
		vector<TrainItem> positive_data_list;
		vector<TrainItem> negative_data_list;
		if (aTargetList.size() != aSVDataList.size())
			throw range_error("ERROR2.20: number of target values: " + stream_cast<string>(aTargetList.size()) + " is different from dataset size:" + stream_cast<string>(aSVDataList.size()));
		for (unsigned i = 0; i < aTargetList.size(); ++i) {
			TrainItem ti;
			ti.mTarget = aTargetList[i];
			ti.mpInstance = aSVDataList[i];
			ti.mpMargin = &oMarginList[i];
			if (i < aDatasetIDList.size()) { //Synthesized instances are appended after the real instances, so the size information of the original id_list marks the start of the syntesized instances
				ti.mInstanceID = aDatasetIDList[i];
			} else
				ti.mInstanceID = -1; //if the instance has been synthesized then it has no correspondent original graph
			if (aTargetList[i] == 1)
				positive_data_list.push_back(ti);
			else if (aTargetList[i] == -1)
				negative_data_list.push_back(ti);
			else
				throw range_error("ERROR2.21: target has to be 1 or -1; cannot be: " + stream_cast<string>(aTargetList[i]));
		}
		//randomly shuffle data
		for (unsigned i = 0; i < positive_data_list.size(); ++i) {
			unsigned j = rand() * positive_data_list.size() / RAND_MAX;
			swap(positive_data_list[i], positive_data_list[j]);
		}
		for (unsigned i = 0; i < negative_data_list.size(); ++i) {
			unsigned j = rand() * negative_data_list.size() / RAND_MAX;
			swap(negative_data_list[i], negative_data_list[j]);
		}

		//over-sample minority class only if there is an imbalance higher than MIN_KFOLD_IMBALANCE and if there is at least one instance for the minority class
		vector<TrainItem> balanced_positive_data_list;
		vector<TrainItem> balanced_negative_data_list;
		const double MIN_KFOLD_IMBALANCE = 1;
		if (p != 0 && p < n / MIN_KFOLD_IMBALANCE) {
			if (!mpParameters->mMinimalOutput)
				cout << "Oversampling positive factor: " << n / (double) p << endl;
			unsigned ratio = n / p;
			unsigned reminder = n % p;
			//duplicate a number of times equal to ratio the datastaset itself
			for (unsigned i = 0; i < ratio; i++)
				balanced_positive_data_list.insert(balanced_positive_data_list.end(), positive_data_list.begin(), positive_data_list.end());
			//add the remainder instances
			for (unsigned i = 0; i < reminder; i++)
				balanced_positive_data_list.push_back(positive_data_list[i]);
			balanced_negative_data_list = negative_data_list;
		} else if (n != 0 && n < p / MIN_KFOLD_IMBALANCE) {
			if (!mpParameters->mMinimalOutput)
				cout << "Oversampling negative factor: " << p / (double) n << endl;
			unsigned ratio = p / n;
			unsigned reminder = p % n;
			for (unsigned i = 0; i < ratio; i++)
				balanced_negative_data_list.insert(balanced_negative_data_list.end(), negative_data_list.begin(), negative_data_list.end());
			for (unsigned i = 0; i < reminder; i++)
				balanced_negative_data_list.push_back(negative_data_list[i]);
			balanced_positive_data_list = positive_data_list;
		} else {
			balanced_positive_data_list = positive_data_list;
			balanced_negative_data_list = negative_data_list;
		}

		//compose dataset by alternating positive and negative examples
		unsigned i;
		for (i = 0; i < balanced_positive_data_list.size(); i++) {
			oDataset.push_back(balanced_positive_data_list[i]);
			if (i < balanced_negative_data_list.size())
				oDataset.push_back(balanced_negative_data_list[i]);
		}
		for (unsigned j = i; j < balanced_negative_data_list.size(); j++)
			oDataset.push_back(balanced_negative_data_list[i]);

		//compute new class ratio
		unsigned bp = 0, bn = 0;
		for (unsigned i = 0; i < oDataset.size(); i++)
			if (oDataset[i].mTarget > 0)
				bp++;
			else
				bn++;
		if (!mpParameters->mMinimalOutput)
			cout << "Rebalanced dataset: " << bp + bn << " (+:" << bp << " -:" << bn << ")" << endl;
	}

	void CoreTrain(vector<TrainItem>& aDataset) {
#ifdef DEBUGON
		VectorClass stat;
#endif

		SVector w_sparsifier;
		SVector w_sparsifier_binarized;
		unsigned num_original_fetures = 0;

		//Iterate epochs times in gradient descent
		ProgressBar pb(1);
		if (!mpParameters->mMinimalOutput) {
			OutputTrainingInfo();
			cout << "Training for " << mpParameters->mEpochs << " epochs." << endl;
		}

		//---------------------------------------------------------------------------------------------------------------
		//iterate for several sparsification iterations
		for (unsigned s = 0; s <= mpParameters->mSparsificationNumIterations; s++) {
			Clear();
			// Shift t in order to have a reasonable initial learning rate. This assumes |x| \approx 1.
			double maxw = 1.0 / sqrt(mpParameters->mLambda);
			double typw = sqrt(maxw);
			double eta0 = typw / max(1.0, dloss(-typw));
			double t = 1 / (eta0 * mpParameters->mLambda);

			if (s != 0 && s == mpParameters->mSparsificationNumIterations) { //for the last iteration (if mSparsificationNumIterations is not just 0 ) use a 0-1 binarized version of the mWSparse
				w_sparsifier_binarized = w_sparsifier;
				w_sparsifier_binarized.binarize();
				if (!mpParameters->mMinimalOutput)
					cout << endl << "Feature filtering step: num features retained: " << w_sparsifier_binarized.sparse_size() << endl;
			}

			//---------------------------------------------------------------------------------------------------------------
			//iterate for several epochs
			for (unsigned e = 0; e < mpParameters->mEpochs; e++) {

				if (e > 0 && s != mpParameters->mSparsificationNumIterations) //after first epoch and excluding last sparsification iteration
					if (mpParameters->mTopologicalRegularizationNumNeighbors != 0) { //if topological regularization is required
						FeatureTopologicalRegularization(1.0 / (mpParameters->mLambda * t));
					}

				if (!mpParameters->mMinimalOutput)
					pb.Count();

				//---------------------------------------------------------------------------------------------------------------
				//iterate over all train instances
				for (unsigned i = 0; i < aDataset.size(); ++i) {
					double eta = 1.0 / (mpParameters->mLambda * t);
					double scale = 1 - eta * mpParameters->mLambda;
					mWScale *= scale;
					if (mWScale < 1e-9) {
						mW.scale(mWScale);
						mWScale = 1;
					}

					//const SVector &x = (*aDataset[i].mpInstance);
					SVector x = (*aDataset[i].mpInstance);

					//iterative sparsification for approximate L0 norm regularization
					if (mpParameters->mSparsificationNumIterations > 0) {
						if (s == mpParameters->mSparsificationNumIterations) { //at last iteration use binarized selection
							assert(w_sparsifier.sparse_size() > 0);
							x = elementwise_product(x, w_sparsifier_binarized);
						} else if (s > 0) {
							assert(w_sparsifier.sparse_size() > 0);
							x = elementwise_product(x, w_sparsifier);
						} else {
							//do not sparsify the very first iteration
						}
					}

					double y = aDataset[i].mTarget;
					double wx = dot(mW, x) * mWScale;
					double margin = (wx + mBias);
					(*(aDataset[i].mpMargin)) = margin;
					double z = y * margin;
#if LOSS < LOGLOSS
					if (z < 1)
#endif
							{
						double etd = eta * dloss(z);
						mW.add(x, etd * y / mWScale);
#if BIAS
						// Slower rate on the bias because it learns at each iteration.
						mBias += etd * y * 0.01;
#endif
					}
					t += 1;
				}
				//---------------------------------------------------------------------------------------------------------------
				//end of iteration over all train instances
				if (s == 0)
					num_original_fetures = mW.sparse_size();
#ifdef DEBUGON
				stat.PushBack((double) mW.sparse_size());
#endif
			}
#ifdef DEBUGON
			cout << endl << "W size statistics: ";
			stat.OutputStatistics(cout);
			cout << endl;
#endif

			if (s != mpParameters->mSparsificationNumIterations) { //skip the last iteration
				//iterative sparsification for approximate L0 norm regularization
				if (s == 0) {
					w_sparsifier = mW;
					w_sparsifier.scale(mWScale);
					//w_sparsifier.normalize();
				} else {
					SVector w_sparsifier_new = mW;
					w_sparsifier_new.scale(mWScale);
					w_sparsifier = elementwise_product(w_sparsifier, w_sparsifier_new, 1e3); //NOTE:set the constant as a hard limit on any elelemnt size
					//w_sparsifier.normalize();
				}
				if (!mpParameters->mMinimalOutput)
					cout << endl << "Iteration: " << s << endl;
				cout << " w_sparsifier norm: " << sqrt(dot(w_sparsifier, w_sparsifier)) << endl;
				cout << " num features retained: " << w_sparsifier.sparse_size() << "/" << num_original_fetures << " (" << (double) (w_sparsifier.sparse_size()) / double(num_original_fetures) << ")" << endl;
			} //end of iterative sparsification

		} //for s
	}

	void FeatureTopologicalRegularization(double aGamma) {
		SVector w_new;
		for (map<unsigned, SVector>::iterator it = mpData->mFeatureCorrelationMatrix.begin(); it != mpData->mFeatureCorrelationMatrix.end(); ++it) {
			unsigned key = it->first;
			SVector& f = it->second;
			double value = dot(f, mW);
			if (value > 1e-9)
				w_new.set(key, value);
		}
		//cout << endl << "w norm before topological regularization:" << sqrt(dot(mW, mW)) * mWScale  << endl; /////////////////
		//cout<<"L*gamma="<<mpParameters->mTopologicalRegularizationRate*aGamma<<endl;//////////////
		mW.combine(1, w_new, -mpParameters->mTopologicalRegularizationRate * aGamma);
		//cout << "w norm after topological regularization:" << sqrt(dot(mW, mW)) * mWScale << endl; /////////////////
	}

	vector<double> Test(vector<unsigned>& aTestSetIDList) {
		vector<double> margin_list;
		for (unsigned i = 0; i < aTestSetIDList.size(); ++i) {
			unsigned gid = aTestSetIDList[i];
			SVector& x = mpData->mVectorList[gid];
			double y = Predict(x);
			margin_list.push_back(y);
		}
		return margin_list;
	}

	inline double Predict(const SVector& x) {
		return dot(mW, x) * mWScale + mBias;
	}

	void OutputTrainingInfo() {
		cout << SEP << endl;
		cout << "Training information" << endl;
		cout << SEP << endl;
		cout << "Lambda: " << mpParameters->mLambda << endl;
		cout << "Epochs: " << mpParameters->mEpochs << endl;
		cout << SEP << endl;
	}

	void OutputModelInfo() {
		cout << SEP << endl;
		cout << "Model information" << endl;
		cout << SEP << endl;
		cout << "W Norm: " << sqrt(dot(mW, mW)) * mWScale << endl;
		cout << "Bias: " << mBias << endl;
		cout << SEP << endl;
	}

};

//------------------------------------------------------------------------------------------------------------------------
class GramMatrixManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
public:
	GramMatrixManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
	}

	void Load() {
		mpData->LoadData();
	}

	void Exec() {
		Load();
		OutputManager();
	}

	void OutputManager() {
		cout << SEP << endl << "Gram matrix phase" << endl << SEP << endl;
		ProgressBar pb;
		pb.Count();
		string output_filename = mpParameters->mInputDataFileName + ".mtx" + mpParameters->mSuffix;
		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.16: Cannot open file:" + output_filename);
		Main(ofs);
		cout << "Gram matrix saved in file " << output_filename << endl;
	}

	void Main(ostream& out) {
		if (!mpParameters->mMinimalOutput){
			cout << "Computing Gram matrix for [" << mpData->mRowIndexList.size()<<" ros x "<<mpData->mColIndexList.size()<<" columns]="<< mpData->mRowIndexList.size() * mpData->mColIndexList.size() << " pairs of instances." << endl;
		}

		double k_sum = 0;
		unsigned counter = 0;
		{
			ProgressBar ppb;

			for (unsigned i = 0; i < mpData->mRowIndexList.size(); ++i) {
				unsigned ii = mpData->mRowIndexList[i];
				for (unsigned j = 0; j < mpData->mColIndexList.size(); ++j) {
					unsigned jj = mpData->mColIndexList[j];
					double k = mpData->ComputeKernel(ii, jj);
					out << k << " ";
					k_sum += k;
					counter++;
				}
				out << endl;
				ppb.Count();
			}
			double avg_k = k_sum / (double) counter;
			cout << endl << "Average kernel value: " << avg_k << endl;
		}
	}
};

//------------------------------------------------------------------------------------------------------------------------
class EmbedManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	double mTau;
	vector<set<unsigned> > mNeighborhoodList;
	vector<set<unsigned> > mNonNeighborhoodList;
public:
	EmbedManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
	}

	void Load() {
		mpData->LoadData();
	}

	void Exec() {
		Load();
		Main();
	}

	void Main() {
		cout << SEP << endl << "Local Multi Dimensional Scaling phase" << endl << SEP << endl;
		vector<FVector> x_list;
		ComputeLocalMultiDimensionalScaling(x_list);

		SaveEmbedding(x_list);
		SaveDistortion(x_list);
		SaveNeighbourhoodList();
	}

	double Norm(const FVector& aX) {
		return sqrt(dot(aX, aX));
	}

	double Distance(const FVector& aX, const FVector& aZ) {
		FVector diff = combine(aX, 1, aZ, -1);
		return Norm(diff);
	}

	FVector Versor(const FVector& aX, const FVector& aZ) {
		FVector diff = combine(aX, 1, aZ, -1);
		diff.scale(Norm(diff));
		return diff;
	}

	void ComputeLocalMultiDimensionalScaling(vector<FVector>& oXList) {
		const unsigned NUM_STEPS_IN_NEIGHBORHOOD_RANGE = 3;
		const double STEP_SIZE_POWER = 1;
		//for random_restart
		vector<FVector> best_x_list;
		double best_distortion = 1;
		double distortion = 1;
		unsigned best_neighborhood_size = 0;
		double best_log_counter = 0;
		cout << "Computing low dimensional layout normalized distortion" << endl;
		unsigned step = (2 * mpParameters->mLMDSNeighborhoodSizeRange / NUM_STEPS_IN_NEIGHBORHOOD_RANGE);
		step = step < 1 ? 1 : step;
		for (unsigned neighborhood_size_modifier = 0; neighborhood_size_modifier <= 2 * mpParameters->mLMDSNeighborhoodSizeRange; neighborhood_size_modifier += step) {
			unsigned effective_neighborhood_size = mpParameters->mLMDSNeighborhoodSize + neighborhood_size_modifier - mpParameters->mLMDSNeighborhoodSizeRange;
			if (effective_neighborhood_size >= 3) {
				InitNeighbourhoodList(effective_neighborhood_size, mpParameters->mLMDSNonNeighborhoodSize);
				for (double log_counter = 0; log_counter <= mpParameters->mLMDSTauExponentRange; log_counter += STEP_SIZE_POWER) {
					double repulsive_force_tau = mTau * mpParameters->mLMDSTau * pow(10, log_counter);
					for (unsigned random_restart_counter = 0; random_restart_counter < mpParameters->mLMDSNumRandomRestarts; random_restart_counter++) {
						cout << "Neighborhood size: " << effective_neighborhood_size << " [" << mpParameters->mLMDSNeighborhoodSize - mpParameters->mLMDSNeighborhoodSizeRange << ".." << mpParameters->mLMDSNeighborhoodSize + mpParameters->mLMDSNeighborhoodSizeRange << "] " << endl << "Repulsive force: "
								<< repulsive_force_tau << endl << "Restart num " << random_restart_counter + 1 << "/" << mpParameters->mLMDSNumRandomRestarts << " (max num iterations: " << mpParameters->mLMDSNumIterations << ")" << endl;
						//initialization phase: random coordinates, computation of long-short distance tradeoff
						vector<FVector> current_x_list;
						LowDimensionalCoordinateInitialization(current_x_list);

						LMDS(current_x_list, repulsive_force_tau);

						distortion = Distortion(current_x_list);
						//keep solution with lowest distortion
						if (distortion < best_distortion) {
							best_x_list = current_x_list;
							best_distortion = distortion;
							best_neighborhood_size = effective_neighborhood_size;
							best_log_counter = log_counter;
							cout << endl << "Saving solution: achieved new low distortion: " << best_distortion << endl;
							SaveEmbedding(best_x_list);
							SaveDistortion(best_x_list);
							SaveNeighbourhoodList();
						} else {
							cout << endl << "Current distortion: " << distortion << endl;
						}
					}
				}
			}
		}
		cout << endl;
		cout << "Best solution found at distortion level: " << best_distortion << " Neighborhood size: " << best_neighborhood_size << " repulsive force multiplicative factor: " << best_log_counter << endl;
		oXList = best_x_list;
	}

	void LMDS(vector<FVector>& current_x_list, double repulsive_force_tau) {
		//iterative minimization loop
		ProgressBar pb;
		for (unsigned itera = 0; itera < mpParameters->mLMDSNumIterations; ++itera) {
			pb.Count();
			for (unsigned i = 0; i < mpData->Size(); ++i) {
				//forall instances in the neighborhood
				for (set<unsigned>::iterator jt = mNeighborhoodList[i].begin(); jt != mNeighborhoodList[i].end(); ++jt) {
					unsigned j = *jt;
					FVector diff_versor = Versor(current_x_list[j], current_x_list[i]);
					double current_distance = Distance(current_x_list[j], current_x_list[i]);
					double desired_distance = (1 - mpData->ComputeKernel(i, j));
					//double stress = sqrt((current_distance - desired_distance) * (current_distance - desired_distance));
					double stress = fabs(current_distance - desired_distance);
					current_x_list[i].add(diff_versor, stress * mpParameters->mLMDSIterationEpsilon);
				}
				//forall instances not in the neighborhood
				for (set<unsigned>::iterator jt = mNonNeighborhoodList[i].begin(); jt != mNonNeighborhoodList[i].end(); ++jt) {
					unsigned j = *jt;
					FVector diff_versor = Versor(current_x_list[i], current_x_list[j]);
					//double long_distance_repulsive_force = mTau * mpParameters->mLMDSTau * (1 - (double) (itera + 1) / (double) mpParameters->mLMDSNumIterations);
					double long_distance_repulsive_force = repulsive_force_tau;
					current_x_list[i].add(diff_versor, long_distance_repulsive_force * mpParameters->mLMDSIterationEpsilon);
				}
			}
		}
	}

	void SaveEmbedding(vector<FVector>& aXList) {
		string output_filename = mpParameters->mInputDataFileName + ".embed" + mpParameters->mSuffix;
		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.45: Cannot open file:" + output_filename);
		for (unsigned t = 0; t < aXList.size(); t++) {
			for (unsigned i = 0; i < mpParameters->mLMDSDimensionality; i++) {
				double val = aXList[t].get(i);
				ofs << val << " ";
			}
			ofs << endl;
		}
		cout << "Embedding saved in file " << output_filename << endl;
	}

	void SaveDistortion(vector<FVector>& aXList) {
		string output_filename = mpParameters->mInputDataFileName + ".distortion" + mpParameters->mSuffix;
		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.56: Cannot open file:" + output_filename);

		vector<set<unsigned> > low_dim_neighbourhood_list;
		MakeNeighbourhoodList(aXList, low_dim_neighbourhood_list);
		assert(aXList.size() == mpData->Size());
		for (unsigned i = 0; i < aXList.size(); ++i) {
			//compute normalized size of the neighborhood intersection
			set<unsigned> intersection;
			set_intersection(mNeighborhoodList[i].begin(), mNeighborhoodList[i].end(), low_dim_neighbourhood_list[i].begin(), low_dim_neighbourhood_list[i].end(), inserter(intersection, intersection.begin()));
			double distortion = 1 - (double) intersection.size() / sqrt((double) mNeighborhoodList[i].size() * (double) low_dim_neighbourhood_list[i].size());
			ofs << distortion << endl;
		}

		cout << "Distortion saved in file " << output_filename << endl;
	}

	void SaveNeighbourhoodList() {
		string output_filename = mpParameters->mInputDataFileName + ".neighbourhood" + mpParameters->mSuffix;
		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.57: Cannot open file:" + output_filename);
		for (unsigned i = 0; i < mpData->Size(); ++i) {
			for (set<unsigned>::iterator it = mNeighborhoodList[i].begin(); it != mNeighborhoodList[i].end(); ++it) {
				ofs << (*it) << " ";
			}
			ofs << endl;
		}
		cout << "Neighbors identity saved in file " << output_filename << endl;
	}

	void InitNeighbourhoodList(unsigned aNeighborhoodSize, unsigned aNonNeighborhoodSize) {
		cout << "Neighborhood indicators computation" << endl;
		mNeighborhoodList.clear();
		mNonNeighborhoodList.clear();
		//for all instances
		vector<double> all_distance_list;
		{
			ProgressBar pb;
			for (unsigned i = 0; i < mpData->Size(); ++i) {
				vector<pair<double, unsigned> > similarity_list(mpData->Size());
				//determine all pairwise similarities
				for (unsigned j = 0; j < mpData->Size(); ++j) {
					double distance = 1 - mpData->ComputeKernel(i, j);
					similarity_list[j] = make_pair(distance, j);
					all_distance_list.push_back(distance);
				}
				//sort and retain the closest k indices
				sort(similarity_list.begin(), similarity_list.end());
				set<unsigned> neighbour_list;
				set<unsigned> non_neighbour_list;
				unsigned effective_size = min(aNeighborhoodSize, (unsigned) similarity_list.size());
				for (unsigned k = 0; k < effective_size; ++k) {
					unsigned neighbour_id = similarity_list[k].second;
					if (neighbour_id != i)
						neighbour_list.insert(neighbour_id);
				}
				unsigned step = (similarity_list.size() - 2 * effective_size) / aNonNeighborhoodSize;
				step = (step == 0) ? 1 : step;
				for (unsigned k = 2 * effective_size; k < similarity_list.size(); k = k + step) {
					unsigned neighbour_id = similarity_list[k].second;
					non_neighbour_list.insert(neighbour_id);
				}
				mNeighborhoodList.push_back(neighbour_list);
				mNonNeighborhoodList.push_back(non_neighbour_list);
				pb.Count();
			}
		}
		//init Tau
		cout << "Computing local-global tradeoff factor: ";
		sort(all_distance_list.begin(), all_distance_list.end());
		unsigned median_index = all_distance_list.size() / 2;
		double median_distance = all_distance_list[median_index];
		mTau = (double) mNeighborhoodList[0].size() / (double) mNonNeighborhoodList[0].size() * median_distance;
		cout << mTau << endl;
	}

	void MakeNeighbourhoodList(vector<FVector>& aXList, vector<set<unsigned> >& oNeighbourhoodList) {
		for (unsigned i = 0; i < aXList.size(); ++i) {
			vector<pair<double, unsigned> > similarity_list(aXList.size());
			//determine all pairwise similarities
			for (unsigned j = 0; j < aXList.size(); ++j) {
				if (i != j) { //NOTE: exclude self
					FVector difference_x = combine(aXList[i], 1, aXList[j], -1);
					double distance = dot(difference_x, difference_x);
					similarity_list[j] = make_pair(distance, j);
				}
			}
			//sort and retain the closest k indices
			sort(similarity_list.begin(), similarity_list.end());
			set<unsigned> neighbour_list;
			unsigned effective_size = min(mpParameters->mLMDSNeighborhoodSize, (unsigned) similarity_list.size());
			for (unsigned k = 0; k < effective_size; ++k) {
				unsigned neighbour_id = similarity_list[k].second;
				neighbour_list.insert(neighbour_id);
			}
			oNeighbourhoodList.push_back(neighbour_list);
		}
	}

	void LowDimensionalCoordinateInitialization(vector<FVector>& oXList) {
		oXList.clear();
		const double span = 2;
		for (unsigned i = 0; i < mpData->Size(); ++i) {
			FVector x;
			for (unsigned j = 0; j < mpParameters->mLMDSDimensionality; ++j) {
				double value = span * random01() - span / 2;
				x.set(j, value);
			}
			oXList.push_back(x);
		}
	}

	double Distortion(vector<FVector>& aXList) {
		unsigned local_continuity = 0;
		vector<set<unsigned> > low_dim_neighbourhood_list;
		MakeNeighbourhoodList(aXList, low_dim_neighbourhood_list);
		assert(aXList.size() == mpData->Size());
		for (unsigned i = 0; i < aXList.size(); ++i) {
			//compute size of the neighborhood intersection
			set<unsigned> intersection;
			set_intersection(mNeighborhoodList[i].begin(), mNeighborhoodList[i].end(), low_dim_neighbourhood_list[i].begin(), low_dim_neighbourhood_list[i].end(), inserter(intersection, intersection.begin()));
			local_continuity += intersection.size();
		}
		double average_local_continuity = ((double) local_continuity / (double) mpParameters->mLMDSNeighborhoodSize) / (double) mpData->Size();
		double average_local_continuity_adjusted_for_chance = average_local_continuity - (double) mpParameters->mLMDSNeighborhoodSize / (double) mpData->Size();
		return 1 - average_local_continuity_adjusted_for_chance;
	}

	double Stress(vector<FVector>& aXList) {
		assert(aXList.size() == mpData->Size());
		double stress = 0;

		//compute average difference of distances
		unsigned counter = 0;
		for (unsigned i = 0; i < aXList.size(); ++i) {
			for (unsigned j = 0; j < aXList.size(); ++j) {
				if (i != j) {
					double desired_distance = 1 - mpData->ComputeKernel(i, j);
					double current_distance = Distance(aXList[i], aXList[j]);
					stress += fabs(current_distance - desired_distance) / (desired_distance);
					//stress += sqrt((current_distance - desired_distance) * (current_distance - desired_distance))/(desired_distance) ;
					counter++;
				}
			}
		}
		double average_stress = stress / (double) counter;
		return average_stress;
	}
};

//------------------------------------------------------------------------------------------------------------------------
class TargetAlignmentManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
public:
	TargetAlignmentManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
	}

	void Load() {
		mpData->LoadData();
	}

	void Exec() {
		Load();
		Main();
	}

	void Main() {
		cout << SEP << endl << "Target alignment phase" << endl << SEP << endl;
		double ta = 0;
		double ka = 0;
		{
			cout << "Computing target alignment for " << mpData->Size() << " instances." << endl;
			ProgressBar ppb;
			for (unsigned i = 0; i < mpData->Size(); ++i) {
				for (unsigned j = 0; j < mpData->Size(); ++j) {
					if (i != j) {
						double k_ij = mpData->ComputeKernel(i, j);
						double t_ij = mpData->mTargetList[i] * mpData->mTargetList[j];
						ta += k_ij * t_ij;
						ka += k_ij * k_ij;
					}
				}
				ppb.Count();
			}
		}
		double target_alignment = ta / sqrt(ka * (mpData->Size() * mpData->Size() - mpData->Size()));
		cout << "Target alignment= " << target_alignment << endl;
	}
};

//------------------------------------------------------------------------------------------------------------------------
class StochasticGradientDescentSupportVectorMachineManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	StochasticGradientDescentSupportVectorMachine mSGDSVM;
public:
	StochasticGradientDescentSupportVectorMachineManager() {
	}

	StochasticGradientDescentSupportVectorMachineManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mSGDSVM.Init(mpParameters, mpData);
	}

	void LoadTarget() {
		mpData->LoadTarget();
	}

	void LoadData() {
		mpData->LoadData();
	}

	void LoadModel() {
		string filename = mpParameters->mModelFileName + mpParameters->mSuffix;
		ifstream ifs;
		ifs.open(filename.c_str());
		if (!ifs)
			throw range_error("ERROR2.23: Cannot open file:" + filename);
		if (!mpParameters->mMinimalOutput)
			cout << endl << "Loading model file: " << filename << endl;
		mSGDSVM.Load(ifs);
		if (!mpParameters->mMinimalOutput) {
			mSGDSVM.OutputModelInfo();
			cout << endl;
		}
	}

	void OutputPerformanceMeasures(ostream& out, const vector<double>& aMarginList, const vector<double>& aTargetList) {
		mSGDSVM.OutputPerformanceMeasures(out, aMarginList, aTargetList);
	}

	double ComputeBalancedFMeasure(const vector<double>& aMarginList, const vector<double>& aTargetList) {
		return mSGDSVM.ComputeBalancedFMeasure(aMarginList, aTargetList);
	}

	double Predict(const SVector& x) {
		return mSGDSVM.Predict(x);
	}

	vector<double> Test(vector<unsigned> aTestSetIDList) {
		vector<unsigned> testset_id_list = aTestSetIDList;
		vector<double> margin_list = mSGDSVM.Test(testset_id_list);
		return margin_list;
	}

	double Test(GraphClass& aG) {
		SVector x;
		mpData->mKernel.GenerateFeatureVector(aG, x);
		return Predict(x);
	}

	vector<double> TestPart(GraphClass& aG) {
		vector<double> margin_list;
		vector<SVector> graph_vertex_vector_list;
		mpData->mKernel.GenerateVertexFeatureVector(aG, graph_vertex_vector_list);
		//for each vertex, compute margin
		unsigned size = mpParameters->mGraphType == "DIRECTED" ? graph_vertex_vector_list.size() / 2 : graph_vertex_vector_list.size();
		for (unsigned vertex_id = 0; vertex_id < size; ++vertex_id) {
			double margin = Predict(graph_vertex_vector_list[vertex_id]);
			if (mpParameters->mGraphType == "DIRECTED")
				margin += Predict(graph_vertex_vector_list[vertex_id + size]);
			margin_list.push_back(margin);
		}
		return margin_list;
	}

	void Train() {
		if (!mpParameters->mMinimalOutput)
			cout << endl << SEP << endl << "Train phase" << endl << SEP << endl;
		ProgressBar pb;
		pb.Count();

		vector<unsigned> train_id_list;
		for (unsigned i = 0; i < mpData->mTargetList.size(); ++i)
			train_id_list.push_back(i);
		Train(mpData->mTargetList, train_id_list);
		SaveModel();

		if (!mpParameters->mMinimalOutput)
			cout << endl << "Train phase completed:";
	}

	void SaveModel(string aLocalSuffix = string()) {
		string filename = mpParameters->mModelFileName + aLocalSuffix + mpParameters->mSuffix;
		ofstream ofs;
		ofs.open(filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.22: Cannot open file:" + filename);
		if (!mpParameters->mMinimalOutput)
			cout << endl << "Saving model file: " << filename << endl;
		mSGDSVM.Save(ofs);
	}

	void Train(vector<double> aTargetList, vector<unsigned> aTrainSetIDList) {
		assert(aTargetList.size() == aTrainSetIDList.size());
		//wrapper for semi-supervised case: self-training
		//assume unsupervised material receives 0 target
		//filter the unsupervised material and put it in separate lists
		//iterate:
		//train on supervised and test on unsupervised
		//replace 0 target with prediction
		vector<double> target_list(aTargetList);
		vector<unsigned> train_supervised_id_list;
		vector<double> train_supervised_target_list;
		vector<unsigned> train_unsupervised_id_list;
		for (unsigned i = 0; i < target_list.size(); ++i) {
			unsigned id = aTrainSetIDList[i];
			double target = target_list[i];
			if (target != 0) {
				train_supervised_id_list.push_back(id);
				train_supervised_target_list.push_back(target);
			} else {
				train_unsupervised_id_list.push_back(id);
			}
		}
		if (train_unsupervised_id_list.size() > 0) {
			//if unsupervised material is present then
			//train on supervised material
			if (!mpParameters->mMinimalOutput) {
				cout << endl << "Semisupervised training on " << target_list.size() << " instances" << endl;
				cout << TAB << "supervised instances: " << train_supervised_id_list.size() << " (" << 100 * train_supervised_id_list.size() / (double) target_list.size() << "%)" << endl;
				cout << TAB << "unsupervised instances: " << train_unsupervised_id_list.size() << " (" << 100 * train_unsupervised_id_list.size() / (double) target_list.size() << "%)" << endl;
			}
			vector<double> margin_list = mSGDSVM.Train(train_supervised_target_list, train_supervised_id_list);

			//repeat for a predefined number of iteration
			for (unsigned iteration = 0; iteration < mpParameters->mSemiSupervisedNumIterations; iteration++) {
				//test on unsupervised material
				if (!mpParameters->mMinimalOutput) {
					cout << endl << TAB << "Iteration " << iteration + 1 << "/" << mpParameters->mSemiSupervisedNumIterations << endl;
					cout << "Testing on unsupervised instances: " << train_unsupervised_id_list.size() << endl;
				}
				vector<double> margin_list = Test(train_unsupervised_id_list);
				//find high and low threshold for margin (i.e. high confidence predictions)
				vector<double> sorted_margin_list;
				vector<double> sorted_positive_margin_list;
				vector<double> sorted_negative_margin_list;
				for (unsigned i = 0; i < margin_list.size(); ++i) {
					sorted_margin_list.push_back(margin_list[i]);
					if (margin_list[i] > 0)
						sorted_positive_margin_list.push_back(margin_list[i]);
					else
						sorted_negative_margin_list.push_back(margin_list[i]);
				}

				unsigned high_threshold_id, low_threshold_id;
				double high_threshold, low_threshold;
				if (!mpParameters->mMinimalOutput)
					cout << "Predicted class distribution:  +:" << sorted_positive_margin_list.size() << " (" << 100 * sorted_positive_margin_list.size() / (double) train_unsupervised_id_list.size() << " %)" << " -:" << sorted_negative_margin_list.size() << " ("
							<< 100 * sorted_negative_margin_list.size() / (double) train_unsupervised_id_list.size() << " %)" << endl;
				if (sorted_positive_margin_list.size() == 0 || sorted_negative_margin_list.size() == 0) {
					if (!mpParameters->mMinimalOutput)
						cout << "Warning: margins are one sided. Proceeding to use margin rank irrespectively of margin sign. Retaining " << mpParameters->mSemiSupervisedThreshold / 2 * 100 << "% of most reliable predictions" << endl;
					sort(sorted_margin_list.begin(), sorted_margin_list.end());
					high_threshold_id = (sorted_margin_list.size() - 1) * (1 - mpParameters->mSemiSupervisedThreshold / 2);
					low_threshold_id = (sorted_margin_list.size() - 1) * (mpParameters->mSemiSupervisedThreshold / 2);
					high_threshold = sorted_margin_list.size() > 0 ? sorted_margin_list[high_threshold_id] : sorted_margin_list[sorted_margin_list.size() - 1];
					low_threshold = sorted_margin_list.size() > 0 ? sorted_margin_list[low_threshold_id] : sorted_margin_list[0];
				} else {
					sort(sorted_positive_margin_list.begin(), sorted_positive_margin_list.end());
					sort(sorted_negative_margin_list.begin(), sorted_negative_margin_list.end());
					high_threshold_id = (sorted_positive_margin_list.size() - 1) * (1 - mpParameters->mSemiSupervisedThreshold);
					low_threshold_id = (sorted_negative_margin_list.size() - 1) * (mpParameters->mSemiSupervisedThreshold);
					high_threshold = sorted_positive_margin_list.size() > 0 ? sorted_positive_margin_list[high_threshold_id] : sorted_negative_margin_list[sorted_negative_margin_list.size() - 1];
					low_threshold = sorted_negative_margin_list.size() > 0 ? sorted_negative_margin_list[low_threshold_id] : sorted_positive_margin_list[0];
				}
				if (!mpParameters->mMinimalOutput)
					cout << "Low score threshold:" << low_threshold << " High score threshold:" << high_threshold << endl;
				//replace 0 target with predicted target only for high confidence predictions
				map<unsigned, double> semi_supervise_augmented_target_map;
				for (unsigned i = 0; i < target_list.size(); ++i) //copy target for supervised instances
					if (target_list[i] != 0) {
						unsigned id = aTrainSetIDList[i];
						semi_supervise_augmented_target_map[id] = target_list[i];
					}
				unsigned counter_p = 0;
				unsigned counter_n = 0;
				for (unsigned i = 0; i < train_unsupervised_id_list.size(); ++i) { //copy prediction for unsupervised instances
					unsigned id = train_unsupervised_id_list[i];
					double margin = margin_list[i];
					double predicted_target;
					bool margin_test = false;
					if (mpParameters->mSemiSupervisedInduceOnlyPositive) {
						margin_test = (margin >= high_threshold);
						predicted_target = 1;
					} else if (mpParameters->mSemiSupervisedInduceOnlyNegative) {
						margin_test = (margin <= low_threshold);
						predicted_target = -1;
					} else {
						margin_test = (margin <= low_threshold || margin >= high_threshold);
						predicted_target = margin <= low_threshold ? -1 : 1;
					}
					if (margin_test == true) {
						assert(semi_supervise_augmented_target_map.count(id) == 0);
						semi_supervise_augmented_target_map[id] = predicted_target;
						if (predicted_target > 0)
							counter_p++;
						else
							counter_n++;
					}
				}
				if (mpParameters->mSemiSupervisedInduceOnlyPositive)
					if (!mpParameters->mMinimalOutput)
						cout << "Adding only predicted positives" << endl;
				if (mpParameters->mSemiSupervisedInduceOnlyNegative)
					if (!mpParameters->mMinimalOutput)
						cout << "Adding only predicted negatives" << endl;
				if (!mpParameters->mMinimalOutput)
					cout << "Added +:" << counter_p << " and -:" << counter_n << " instances from unsupervised set to training set of size " << train_supervised_id_list.size() << endl;
				//compose indices vectors for training instances and target
				vector<unsigned> train_semi_supervise_augmented_id_list;
				vector<double> train_semi_supervise_augmented_target_list;
				for (map<unsigned, double>::iterator it = semi_supervise_augmented_target_map.begin(); it != semi_supervise_augmented_target_map.end(); ++it) {
					unsigned id = it->first;
					double target = it->second;
					train_semi_supervise_augmented_id_list.push_back(id);
					train_semi_supervise_augmented_target_list.push_back(target);
				}
				//retrain
				margin_list = mSGDSVM.Train(train_semi_supervise_augmented_target_list, train_semi_supervise_augmented_id_list);
			}

		} else { //if no unsupervised material is present then train directly
			vector<double> margin_list = mSGDSVM.Train(target_list, aTrainSetIDList);
		}
	}
};

//------------------------------------------------------------------------------------------------------------------------
class CrossValidationManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	StochasticGradientDescentSupportVectorMachineManager mSGDSVMManager;
public:
	CrossValidationManager() {
	}

	CrossValidationManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mSGDSVMManager.Init(apParameters, apData);
	}

	void Load() {
		mSGDSVMManager.LoadTarget();
		mSGDSVMManager.LoadData();
	}

	void Exec() {
		Load();
		Main();
	}

	void Main() {
		map<unsigned, vector<double> > results;
		CrossValidation(results);
	}

	double GetBalancedFMeasure() {
		map<unsigned, vector<double> > results;
		return CrossValidation(results);
	}

	double CrossValidation(map<unsigned, vector<double> >& oTestResultMap) {
		oTestResultMap.clear();
		srand(mpParameters->mRandomSeed);
		ProgressBar pbt;
		pbt.Count();
		vector<pair<double, double> > prediction_list;
		vector<pair<double, double> > margin_list;

		//randomly shuffle indices
		unsigned size = mpData->Size();
		vector<unsigned> data_id_list;
		for (unsigned i = 0; i < size; ++i)
			data_id_list.push_back(i);
		for (unsigned i = 0; i < size; ++i) {
			unsigned j = rand() * size / RAND_MAX;
			swap(data_id_list[i], data_id_list[j]);
		}

		//loop to build train-test split in the cross validation way
		for (unsigned f = 0; f < mpParameters->mCrossValidationNumFolds; f++) {
			ProgressBar pbf;
			pbf.Count();
			if (!mpParameters->mMinimalOutput)
				cout << SEP << endl << TAB << TAB << "Fold: " << f + 1 << " of " << mpParameters->mCrossValidationNumFolds << endl;
			vector<unsigned> train_id_list;
			vector<unsigned> test_id_list;
			map<unsigned, unsigned> class_counter_map;
			for (unsigned i = 0; i < size; ++i) {
				unsigned id = data_id_list[i];
				double target = mpData->mTargetList[id];
				if (class_counter_map.count(target) == 0)
					class_counter_map[target] = 1;
				else
					class_counter_map[target]++;
				if (target != 0 && class_counter_map[target] % mpParameters->mCrossValidationNumFolds == f) //NOTE: exclude un supervised material from test element list
					test_id_list.push_back(id);
				else
					train_id_list.push_back(id);
			}

			//extract target list for training
			vector<double> train_target_list;
			for (unsigned i = 0; i < train_id_list.size(); i++) {
				unsigned id = train_id_list[i];
				train_target_list.push_back(mpData->mTargetList[id]);
			}
			//extract target list for testing
			vector<double> test_target_list;
			for (unsigned i = 0; i < test_id_list.size(); i++) {
				unsigned id = test_id_list[i];
				test_target_list.push_back(mpData->mTargetList[id]);
			}

			//perform training
			mSGDSVMManager.Train(train_target_list, train_id_list);
			string model_filename_prefix = "_" + stream_cast<string>(f + 1);
			mSGDSVMManager.SaveModel(model_filename_prefix);

			//perform testing
			vector<double> fold_margin_list = mSGDSVMManager.Test(test_id_list);
			//add to test_result_map
			assert(fold_margin_list.size() == test_id_list.size());
			for (unsigned i = 0; i < fold_margin_list.size(); i++) {
				unsigned id = test_id_list[i];
				//pack all the result fields sequentially in a vector
				double target = test_target_list[i];
				double margin = fold_margin_list[i];
				double prediction = margin > 0 ? 1 : -1;
				vector<double> res;
				res.push_back(target);
				res.push_back(prediction);
				res.push_back(margin);
				//memoize the result vector with the test id
				oTestResultMap[id] = res;
			}
			if (!mpParameters->mMinimalOutput)
				cout << "Fold phase concluded in:" << endl;
		}

		vector<double> cv_prediction_list;
		vector<double> cv_target_list;
		string ofs_name = mpParameters->mInputDataFileName + ".cv_predictions" + mpParameters->mSuffix;
		ofstream ofs(ofs_name.c_str());
		//for all test ids read in order
		for (map<unsigned, vector<double> >::iterator it = oTestResultMap.begin(); it != oTestResultMap.end(); ++it) {
			unsigned id = it->first;
			//unpack the result fields from the result vector memoized with the test id
			double target = it->second[0];
			double prediction = it->second[1];
			double margin = it->second[2];
			ofs << id << " " << target << " " << prediction << " " << margin << endl;
			cv_prediction_list.push_back(prediction);
			cv_target_list.push_back(target);
		}
		if (!mpParameters->mMinimalOutput) {
			cout << SEP << endl << "Performance on data set in cross validation:" << endl;
			mSGDSVMManager.OutputPerformanceMeasures(cout, cv_prediction_list, cv_target_list);
			cout << endl << "Instance id, true target, prediction and margin saved in file: " << ofs_name << endl;
			cout << "Crossvalidation concluded in:" << endl;
		}
		double bfmeasure = mSGDSVMManager.ComputeBalancedFMeasure(cv_prediction_list, cv_target_list);
		return bfmeasure;
	}

};

//------------------------------------------------------------------------------------------------------------------------
class LearningCurveManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	StochasticGradientDescentSupportVectorMachineManager mSGDSVMManager;

public:
	LearningCurveManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mSGDSVMManager.Init(apParameters, apData);
	}

	void Load() {
		mSGDSVMManager.LoadTarget();
		mSGDSVMManager.LoadData();
	}

	void Exec() {
		Load();
		Main();
	}

	void Main() {
		LearningCurve();
	}

	void LearningCurve() {
		srand(mpParameters->mRandomSeed);
		cout << SEP << endl;
		cout << "Computing learning curve with " << mpParameters->mNumPoints << " folds." << endl;
		ProgressBar pbt;
		pbt.Count();
		vector<pair<double, double> > prediction_list;
		vector<pair<double, double> > margin_list;

		unsigned size = mpData->Size();
		//randomly shuffle indices
		vector<unsigned> data_id_list;
		for (unsigned i = 0; i < size; ++i)
			data_id_list.push_back(i);
		for (unsigned i = 0; i < size; ++i) {
			unsigned j = (double) rand() / RAND_MAX * size;
			swap(data_id_list[i], data_id_list[j]);
		}

		//build train-test split in the learning curve way
		vector<unsigned> test_id_list;
		//test data is the first fold
		for (unsigned i = 0; i < size / mpParameters->mNumPoints; ++i) {
			unsigned id = data_id_list[i];
			test_id_list.push_back(id);
		}
		//sort the indices in order to guarantee sequential file access
		sort(test_id_list.begin(), test_id_list.end());
		//extract target list for testing
		vector<double> test_target_list;
		for (unsigned i = 0; i < test_id_list.size(); i++) {
			unsigned id = test_id_list[i];
			test_target_list.push_back(mpData->mTargetList[id]);
		}
		//training data is built incrementally adding 1/mpParameters->mLearningCurveNumPoints * size instances
		for (unsigned f = 1; f < mpParameters->mNumPoints; f++) { //NOTE: start from 1 as the first fold is used for the test data
			ProgressBar pbf;
			pbf.Count();
			cout << SEP << endl;
			cout << TAB << TAB << "Fold: " << f << " of " << mpParameters->mNumPoints - 1 << endl;
			//generate the training set
			vector<unsigned> train_id_list;
			for (unsigned i = size / mpParameters->mNumPoints; i < size * (f + 1) / mpParameters->mNumPoints; ++i) {
				unsigned id = data_id_list[i];
				train_id_list.push_back(id);
			}
			//sort the indices in order to guarantee sequential file access
			sort(train_id_list.begin(), train_id_list.end());
			//extract target list for training
			vector<double> train_target_list;
			for (unsigned i = 0; i < train_id_list.size(); i++) {
				unsigned id = train_id_list[i];
				train_target_list.push_back(mpData->mTargetList[id]);
			}

			//perform training
			mSGDSVMManager.Train(train_target_list, train_id_list);

			//compute predictions on test
			{
				vector<double> fold_margin_list = mSGDSVMManager.Test(test_id_list);
				assert(fold_margin_list.size() == test_id_list.size());
				cout << SEP << endl << "Performance on test set:" << endl;
				mSGDSVMManager.OutputPerformanceMeasures(cout, fold_margin_list, test_target_list);

				//save results to file
				string ofs_name = mpParameters->mInputDataFileName + ".lc_predictions_test_fold_" + stream_cast<string>(f) + mpParameters->mSuffix;
				ofstream ofs(ofs_name.c_str());
				for (unsigned i = 0; i < fold_margin_list.size(); i++) {
					unsigned id = test_id_list[i];
					double target = test_target_list[i];
					double margin = fold_margin_list[i];
					double prediction = margin > 0 ? 1 : -1;
					ofs << id << " " << target << " " << prediction << " " << margin << endl;
				}
				ofs.close();
				cout << endl << "Instance id, true target, prediction and margin saved in file: " << ofs_name << endl;
			}
			//compute predictions on train
			{
				vector<double> fold_margin_list = mSGDSVMManager.Test(train_id_list);
				assert(fold_margin_list.size() == train_id_list.size());
				cout << SEP << endl << "Performance on train set:" << endl;
				mSGDSVMManager.OutputPerformanceMeasures(cout, fold_margin_list, train_target_list);

				//save results to file
				string ofs_name = mpParameters->mInputDataFileName + ".lc_predictions_train_fold_" + stream_cast<string>(f) + mpParameters->mSuffix;
				ofstream ofs(ofs_name.c_str());
				for (unsigned i = 0; i < fold_margin_list.size(); i++) {
					unsigned id = train_id_list[i];
					double target = train_target_list[i];
					double margin = fold_margin_list[i];
					double prediction = margin > 0 ? 1 : -1;
					ofs << id << " " << target << " " << prediction << " " << margin << endl;
				}
				ofs.close();
				cout << endl << "Instance id, true target, prediction and margin saved in file: " << ofs_name << endl;
			}

			cout << "Fold phase concluded in:" << endl;
		}
		cout << "Learning curve concluded in:" << endl;
	}
};

//------------------------------------------------------------------------------------------------------------------------
class BiasVarianceDecompositionManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	CrossValidationManager mCrossValidationManager;
public:
	BiasVarianceDecompositionManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mCrossValidationManager.Init(apParameters, apData);
	}

	void Load() {
		mCrossValidationManager.Load();
	}

	void Exec() {
		Load();
		Main();
	}

	void Main() {
		BiasVarianceComputation();
	}

	void BiasVarianceComputation() {
		string output_filename = mpParameters->mInputDataFileName + ".conf" + mpParameters->mSuffix;
		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.45: Cannot open file:" + output_filename);

		ProgressBar pbt;
		pbt.Count();
		vector<map<unsigned, vector<double> > > multi_result_list;
		for (unsigned i = 0; i < mpParameters->mNumPoints; ++i) {
			mpParameters->mRandomSeed += i;
			map<unsigned, vector<double> > result_list;
			mCrossValidationManager.CrossValidation(result_list);
			multi_result_list.push_back(result_list);
		}

		unsigned size = mpData->Size();
		vector<double> accuracy_list(size);
		vector<VectorClass> margin_list(size);
		for (unsigned r = 0; r < mpParameters->mNumPoints; ++r) {
			for (unsigned i = 0; i < size; ++i) {
				double target = multi_result_list[r][i][0];
				double prediction = multi_result_list[r][i][1];
				double margin = multi_result_list[r][i][2];
				margin_list[i].PushBack(margin);
				if (prediction == target)
					accuracy_list[i] += 1 / (double) mpParameters->mNumPoints;
			}
		}

		//write results to file
		for (unsigned i = 0; i < size; ++i) {
			double target = multi_result_list[0][i][0];
			double mean_margin = margin_list[i].Mean();
			double std_margin = margin_list[i].StandardDeviation();
			ofs << target << " " << accuracy_list[i] << " " << mean_margin << " " << std_margin << endl;
		}
		cout << "Target, average accuracy, average margin, margin standard deviation saved in file " << output_filename << endl;
	}

};

//------------------------------------------------------------------------------------------------------------------------
class ParametersOptimizationManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	CrossValidationManager mCrossValidationManager;

	double mCurrentBFMeasure;
	double mBestBFMeasure;

	double mLambdaLimit;
	unsigned mEpochsLimit;
	unsigned mRadiusLimit;
	unsigned mDistanceLimit;
	unsigned mTopologicalRegularizationNumNeighborsLimit;
	double mTopologicalRegularizationRateLimit;
	unsigned mSparsificationNumIterationsLimit;
	double mTreeLambdaLimit;
	double mRadiusTwoLimit;

public:
	ParametersOptimizationManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mCrossValidationManager.Init(apParameters, apData);

		mCurrentBFMeasure = 0;
		mBestBFMeasure = 0;
		mLambdaLimit = mpParameters->mLambda;
		mEpochsLimit = mpParameters->mEpochs;
		mRadiusLimit = mpParameters->mRadius;
		mDistanceLimit = mpParameters->mDistance;
		mTopologicalRegularizationNumNeighborsLimit = mpParameters->mTopologicalRegularizationNumNeighbors;
		mTopologicalRegularizationRateLimit = mpParameters->mTopologicalRegularizationRate;
		mSparsificationNumIterationsLimit = mpParameters->mSparsificationNumIterations;
		mTreeLambdaLimit = mpParameters->mTreeLambda;
		mRadiusTwoLimit = mpParameters->mRadiusTwo;


	}

	void Exec() {
		mpData->LoadTarget();
		OutputManager();
	}

	void OutputManager() {
		string output_filename = mpParameters->mInputDataFileName + ".opt_param" + mpParameters->mSuffix;
		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.45: Cannot open file:" + output_filename);

		ParametersOptimization(ofs);

		cout << "Best parameter configuration obtains a balanced F-measure= " << mBestBFMeasure << endl;
		cout << "For best predictive performance use the following parameter setting:" << endl;
		OutputParameters(cout);
		OutputParameters(ofs);
		cout << "Optimal parameters configuration saved in file " << output_filename << endl;
	}

	void ParametersOptimization(ostream& ofs) {
		ProgressBar pbt;
		pbt.Count();

		SetDefaultParameters();
		cout << "Initial parameters configuration: ";
		OutputParameters(cout);
		mpData->LoadData();

		cout << "Line search parameters optimization: " << mpParameters->mNumLineSearchIterations << " iterations." << endl;
		for (unsigned line_search_iteration = 0; line_search_iteration < mpParameters->mNumLineSearchIterations; line_search_iteration++) {
			cout << "Iteration " << line_search_iteration + 1 << "/" << mpParameters->mNumLineSearchIterations << endl;
			OptimizeLambda();
			OptimizeEpochs();
			OptimizeDistance();
			OptimizeRadius();
			OptimizeTopologicalRegularizationNumNeighbors();
			OptimizeTopologicalRegularizationRate();
			OptimizeSparsificationNumIterations();
			if (mpParameters->mKernelType == "DDK" || mpParameters->mKernelType == "NSDDK" || mpParameters->mKernelType == "ANSDDK"){
				OptimizeTreeLambda();
			}
			if (mpParameters->mKernelType == "NSDDK" || mpParameters->mKernelType == "ANSDDK"){
				OptimizeRadiusTwo();
			}

		}
	}

	//Required for DDK kernel family
	void OptimizeTreeLambda() {
			mBestBFMeasure = 0;
			double lambda_best = mpParameters->mTreeLambda;
			for (double l = 0.2; l <= mTreeLambdaLimit; l= l+0.2) {
				mpParameters->mTreeLambda = l;
				mpData->LoadData();
				mCurrentBFMeasure = mCrossValidationManager.GetBalancedFMeasure();
				if (mCurrentBFMeasure > mBestBFMeasure) {
					mBestBFMeasure = mCurrentBFMeasure;
					lambda_best = l;
					cout << "*";
				}
				cout << "tree_lambda:" << l << " bf:" << mCurrentBFMeasure << " " << endl;
			}
			mpParameters->mTreeLambda = lambda_best;
			//mpData->LoadData();
			cout << "bFmeasure:" << mBestBFMeasure << " current best parameters configuration:";
			OutputParameters(cout);
		}
	void OptimizeRadiusTwo() {
		mBestBFMeasure = 0;
		unsigned r_best = mpParameters->mRadiusTwo;
		for (unsigned r = 0; r <= mRadiusTwoLimit; r++) {
			mpParameters->mRadiusTwo = r;
			mpData->LoadData();
			mCurrentBFMeasure = mCrossValidationManager.GetBalancedFMeasure();
			if (mCurrentBFMeasure > mBestBFMeasure) {
				mBestBFMeasure = mCurrentBFMeasure;
				r_best = r;
				cout << "*";
			}
			cout << "radius_two:" << r << " bf:" << mCurrentBFMeasure << " " << endl;
		}
		mpParameters->mRadiusTwo = r_best;
		cout << "bFmeasure:" << mBestBFMeasure << " current best parameters configuration:";
		OutputParameters(cout);
	}
	//-----------
	void SetDefaultParameters() {
		//default values
		mpParameters->mLambda = 1e-4;
		mpParameters->mEpochs = 5;
		mpParameters->mRadius = 1;
		mpParameters->mDistance = 2;
		mpParameters->mTopologicalRegularizationNumNeighbors = 0;
		mpParameters->mTopologicalRegularizationRate = 0.001;
		mpParameters->mSparsificationNumIterations = 0;
		mpParameters->mRadiusTwo = 1;
		mpParameters->mTreeLambda = 1.2;

	}

	void OptimizeLambda() {
		mBestBFMeasure = 0;
		double lambda_best = mpParameters->mLambda;
		const double LAMBDA_UPPER_BOUND = 0.01;
		//double lambda_step = exp((log(LAMBDA_UPPER_BOUND) - log(lambda_limit)) / (2 * mpParameters->mLearningCurveNumPoints));
		double lambda_step = 10;
		for (double lambda = mLambdaLimit; lambda <= LAMBDA_UPPER_BOUND; lambda *= lambda_step) {
			mpParameters->mLambda = lambda;
			mCurrentBFMeasure = mCrossValidationManager.GetBalancedFMeasure();
			if (mCurrentBFMeasure > mBestBFMeasure) {
				mBestBFMeasure = mCurrentBFMeasure;
				lambda_best = lambda;
				cout << "*";
			}
			cout << "l:" << lambda << " bf:" << mCurrentBFMeasure << " " << endl;
		}
		mpParameters->mLambda = lambda_best;
		cout << "bFmeasure:" << mBestBFMeasure << " current best parameters configuration:";
		OutputParameters(cout);
	}

	void OptimizeEpochs() {
		mBestBFMeasure = 0;
		unsigned epochs_best = mpParameters->mEpochs;
		unsigned epochs_step = (double) mEpochsLimit / (double) mpParameters->mNumPoints;
		epochs_step = epochs_step == 0 ? 1 : epochs_step;
		for (unsigned epochs = epochs_step; epochs <= mEpochsLimit; epochs += epochs_step) {
			mpParameters->mEpochs = epochs;
			mCurrentBFMeasure = mCrossValidationManager.GetBalancedFMeasure();
			if (mCurrentBFMeasure > mBestBFMeasure) {
				mBestBFMeasure = mCurrentBFMeasure;
				epochs_best = epochs;
				cout << "*";
			}
			cout << "e:" << epochs << " bf:" << mCurrentBFMeasure << " " << endl;
		}
		mpParameters->mEpochs = epochs_best;
		cout << "bFmeasure:" << mBestBFMeasure << " current best parameters configuration:";
		OutputParameters(cout);
	}

	void OptimizeDistance() {
		mBestBFMeasure = 0;
		unsigned d_best = mpParameters->mDistance;
		for (unsigned d = 0; d <= mDistanceLimit; d++) {
			mpParameters->mDistance = d;
			mpData->LoadData();
			mCurrentBFMeasure = mCrossValidationManager.GetBalancedFMeasure();
			if (mCurrentBFMeasure > mBestBFMeasure) {
				mBestBFMeasure = mCurrentBFMeasure;
				d_best = d;
				cout << "*";
			}
			cout << "d:" << d << " bf:" << mCurrentBFMeasure << " " << endl;
		}
		mpParameters->mDistance = d_best;
		cout << "bFmeasure:" << mBestBFMeasure << " current best parameters configuration:";
		OutputParameters(cout);
	}

	void OptimizeRadius() {
		mBestBFMeasure = 0;
		unsigned r_best = mpParameters->mRadius;
		for (unsigned r = 0; r <= mRadiusLimit; r++) {
			mpParameters->mRadius = r;
			mpData->LoadData();
			mCurrentBFMeasure = mCrossValidationManager.GetBalancedFMeasure();
			if (mCurrentBFMeasure > mBestBFMeasure) {
				mBestBFMeasure = mCurrentBFMeasure;
				r_best = r;
				cout << "*";
			}
			cout << "r:" << r << " bf:" << mCurrentBFMeasure << " " << endl;
		}
		mpParameters->mRadius = r_best;
		cout << "bFmeasure:" << mBestBFMeasure << " current best parameters configuration:";
		OutputParameters(cout);
	}

	void OptimizeTopologicalRegularizationNumNeighbors() {
		mBestBFMeasure = 0;
		unsigned topological_regularization_num_neighbors_best = mpParameters->mTopologicalRegularizationNumNeighbors;
		unsigned topological_regularization_num_neighbors_step = (double) mTopologicalRegularizationNumNeighborsLimit / (double) mpParameters->mNumPoints;
		topological_regularization_num_neighbors_step = topological_regularization_num_neighbors_step == 0 ? 1 : topological_regularization_num_neighbors_step;
		for (unsigned topological_regularization_num_neighbors = 0; topological_regularization_num_neighbors <= mTopologicalRegularizationNumNeighborsLimit; topological_regularization_num_neighbors += topological_regularization_num_neighbors_step) {
			mpParameters->mTopologicalRegularizationNumNeighbors = topological_regularization_num_neighbors;
			mpData->LoadData();
			mCurrentBFMeasure = mCrossValidationManager.GetBalancedFMeasure();
			if (mCurrentBFMeasure > mBestBFMeasure) {
				mBestBFMeasure = mCurrentBFMeasure;
				topological_regularization_num_neighbors_best = topological_regularization_num_neighbors;
				cout << "*";
			}
			cout << "C:" << topological_regularization_num_neighbors << " bf:" << mCurrentBFMeasure << " " << endl;
		}
		mpParameters->mTopologicalRegularizationNumNeighbors = topological_regularization_num_neighbors_best;
		cout << "bFmeasure:" << mBestBFMeasure << " current best parameters configuration:";
		OutputParameters(cout);
	}

	void OptimizeTopologicalRegularizationRate() {
		mBestBFMeasure = 0;
		double topological_regularization_rate_best = mpParameters->mTopologicalRegularizationRate;
		const double TOPOLOGICAL_REGULARIZATION_RATE_UPPER_BOUND = .001;
		//double topological_regularization_rate_step = (double) topological_regularization_rate_limit / (double) mpParameters->mLearningCurveNumPoints;
		double topological_regularization_rate_step = 10;
		for (double topological_regularization_rate = mTopologicalRegularizationRateLimit; topological_regularization_rate <= TOPOLOGICAL_REGULARIZATION_RATE_UPPER_BOUND; topological_regularization_rate *= topological_regularization_rate_step) {
			mpParameters->mTopologicalRegularizationRate = topological_regularization_rate;
			mCurrentBFMeasure = mCrossValidationManager.GetBalancedFMeasure();
			if (mCurrentBFMeasure > mBestBFMeasure) {
				mBestBFMeasure = mCurrentBFMeasure;
				topological_regularization_rate_best = topological_regularization_rate;
				cout << "*";
			}
			cout << "L:" << topological_regularization_rate << " bf:" << mCurrentBFMeasure << " " << endl;
		}
		mpParameters->mTopologicalRegularizationRate = topological_regularization_rate_best;
		cout << "bFmeasure:" << mBestBFMeasure << " current best parameters configuration:";
		OutputParameters(cout);
	}

	void OptimizeSparsificationNumIterations() {
		mBestBFMeasure = 0;
		double sparsification_num_iterations_best = mpParameters->mSparsificationNumIterations;
		unsigned sparsification_num_iterations_step = (double) mSparsificationNumIterationsLimit / (double) mpParameters->mNumPoints;
		sparsification_num_iterations_step = sparsification_num_iterations_step == 0 ? 1 : sparsification_num_iterations_step;
		for (unsigned sparsification_num_iterations = 0; sparsification_num_iterations <= mSparsificationNumIterationsLimit; sparsification_num_iterations += sparsification_num_iterations_step) {
			mpParameters->mSparsificationNumIterations = sparsification_num_iterations;
			mCurrentBFMeasure = mCrossValidationManager.GetBalancedFMeasure();
			if (mCurrentBFMeasure > mBestBFMeasure) {
				mBestBFMeasure = mCurrentBFMeasure;
				sparsification_num_iterations_best = sparsification_num_iterations;
				cout << "*";
			}
			cout << "O:" << sparsification_num_iterations << " bf:" << mCurrentBFMeasure << " " << endl;
		}
		mpParameters->mSparsificationNumIterations = sparsification_num_iterations_best;
		cout << "bFmeasure:" << mBestBFMeasure << " current best parameters configuration:";
		OutputParameters(cout);
	}

	void OutputParameters(ostream& out) {
		out << " -r " << mpParameters->mRadius;
		out << " -d " << mpParameters->mDistance;
		out << " -e " << mpParameters->mEpochs;
		out << " -l " << mpParameters->mLambda;
		out << " -O " << mpParameters->mSparsificationNumIterations;
		out << " -L " << mpParameters->mTopologicalRegularizationRate;
		out << " -C " << mpParameters->mTopologicalRegularizationNumNeighbors;
		if (mpParameters->mKernelType == "DDK" || mpParameters->mKernelType == "NSDDK" || mpParameters->mKernelType == "ANSDDK")
			out << " --tree_lambda " << mpParameters->mTreeLambda;
		if (mpParameters->mKernelType == "DDK" || mpParameters->mKernelType == "NSDDK" || mpParameters->mKernelType == "ANSDDK")
			out << " --radius_two " << mpParameters->mRadiusTwo;


		out << endl;
	}
};

//------------------------------------------------------------------------------------------------------------------------
class TestManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	StochasticGradientDescentSupportVectorMachineManager mSGDSVMManager;
public:
	TestManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mpData->mKernel.ParametersSetup();
		mSGDSVMManager.Init(apParameters, apData);
	}

	void Load() {
		mSGDSVMManager.LoadModel();
	}

	void Exec() {
		Load();
		InputOutputManager();
	}

	void InputOutputManager() {
		//output
		//compose output filename
		string output_filename = mpParameters->mInputDataFileName;
		output_filename += ".prediction";
		output_filename += mpParameters->mSuffix;

		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.16: Cannot open file: " + output_filename);

		//input
		igzstream fin;
		fin.open(mpParameters->mInputDataFileName.c_str());
		if (!fin)
			throw range_error("ERROR2.11: Cannot open file: " + mpParameters->mInputDataFileName);
		//perform online action
		if (!mpParameters->mMinimalOutput)
			cout << "Processing I/O file: " << mpParameters->mInputDataFileName << endl;
		Main(fin, ofs);
		if (!mpParameters->mMinimalOutput)
			cout << "Result saved in file " + output_filename << endl;
	}

	void Main(istream& fin, ostream& ofs) {
		ProgressBar pb;
		unsigned instance_counter = 0;
		bool valid_input = true;
		while (!fin.eof() && valid_input) {

			switch (mpParameters->mFileTypeCode) {
			case GRAPH:
#ifdef USEOBABEL
			case MOLECULAR_GRAPH:
#endif
				case SEQUENCE: {
				GraphClass g;
				mpData->SetGraphFromFile(fin, g);
				if (!g.IsEmpty()) {
					SVector x;
					mpData->mKernel.GenerateFeatureVector(g, x);
					double margin = mSGDSVMManager.Predict(x);
					int prediction = margin > 0 ? 1 : -1;
					ofs << prediction << " " << margin << endl;
				} else
					valid_input = false;
				if (!mpParameters->mMinimalOutput)
					pb.Count();
				instance_counter++;
			}
				break;
			case SPARSE_VECTOR: {
				SVector x;
				if (mpParameters->mBinaryFormat)
					mpData->SetVectorFromSparseVectorBinaryFile(fin, x);
				else
					mpData->SetVectorFromSparseVectorAsciiFile(fin, x);
				if (x.size() > 0) {
					double margin = mSGDSVMManager.Predict(x);
					int prediction = margin > 0 ? 1 : -1;
					ofs << prediction << " " << margin << endl;
				}
				if (!mpParameters->mMinimalOutput)
					pb.Count();
				instance_counter++;
			}
				break;
			default:
				throw range_error("ERROR2.45: file type not recognized: " + mpParameters->mFileType);
			}
		}
	}
};

//------------------------------------------------------------------------------------------------------------------------
class TestPartManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	StochasticGradientDescentSupportVectorMachineManager mSGDSVMManager;
public:
	TestPartManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mpData->mKernel.ParametersSetup();
		mSGDSVMManager.Init(apParameters, apData);
	}

	void Load() {
		mSGDSVMManager.LoadModel();
	}

	void Exec() {
		Load();
		InputOutputManager();
	}

	void InputOutputManager() {
		//output
		//compose output filename
		string output_filename = mpParameters->mInputDataFileName;
		output_filename += ".prediction_part";
		output_filename += mpParameters->mSuffix;

		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.16: Cannot open file: " + output_filename);

		//input
		igzstream fin;
		fin.open(mpParameters->mInputDataFileName.c_str());
		if (!fin)
			throw range_error("ERROR2.11: Cannot open file: " + mpParameters->mInputDataFileName);
		//perform online action
		if (!mpParameters->mMinimalOutput)
			cout << "Processing I/O file: " << mpParameters->mInputDataFileName << endl;
		Main(fin, ofs);
		if (!mpParameters->mMinimalOutput)
			cout << "Result saved in file " + output_filename << endl;
	}

	void Main(istream& fin, ostream& ofs) {
		ProgressBar pb;

		unsigned instance_counter = 0;
		bool valid_input = true;
		while (!fin.eof() && valid_input) {

			switch (mpParameters->mFileTypeCode) {
			case GRAPH:
#ifdef USEOBABEL
			case MOLECULAR_GRAPH:
#endif
				case SEQUENCE: {
				GraphClass g;
				mpData->SetGraphFromFile(fin, g);
				if (!g.IsEmpty()) {
					vector<SVector> graph_vertex_vector_list;
					mpData->mKernel.GenerateVertexFeatureVector(g, graph_vertex_vector_list);
					//for each vertex, compute margin
					unsigned size = mpParameters->mGraphType == "DIRECTED" ? graph_vertex_vector_list.size() / 2 : graph_vertex_vector_list.size();
					for (unsigned vertex_id = 0; vertex_id < size; ++vertex_id) {
						double margin = mSGDSVMManager.Predict(graph_vertex_vector_list[vertex_id]);
						if (mpParameters->mGraphType == "DIRECTED")
							margin += mSGDSVMManager.Predict(graph_vertex_vector_list[vertex_id + size]);
						ofs << instance_counter << " " << vertex_id << " " << margin << endl;
					}
					if (!mpParameters->mMinimalOutput)
						pb.Count();
					instance_counter++;
				} else
					valid_input = false;
			}
				break;
			case SPARSE_VECTOR:
				throw range_error("ERROR2.4512: Cannot process directly sparse vecotor representations for TEST_PART");
				break;
			default:
				throw range_error("ERROR2.45: file type not recognized: " + mpParameters->mFileType);
			}
		}
	}
};

//------------------------------------------------------------------------------------------------------------------------
class FeatureManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
public:
	FeatureManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mpData->mKernel.ParametersSetup();
	}

	void Exec() {
		InputOutputManager();
	}

	void InputOutputManager() {
		//output
		//compose output filename
		string output_filename = mpParameters->mInputDataFileName;
		output_filename += ".feature";
		output_filename += mpParameters->mSuffix;

		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.16: Cannot open file: " + output_filename);

		//input
		igzstream fin;
		fin.open(mpParameters->mInputDataFileName.c_str());
		if (!fin)
			throw range_error("ERROR2.11: Cannot open file: " + mpParameters->mInputDataFileName);
		//perform online action
		if (!mpParameters->mMinimalOutput)
			cout << "Processing I/O file: " << mpParameters->mInputDataFileName << endl;
		Main(fin, ofs);
		if (!mpParameters->mMinimalOutput)
			cout << "Result saved in file " + output_filename << endl;
	}

	void Main(istream& fin, ostream& ofs) {
		ProgressBar pb;
		unsigned instance_counter = 0;
		bool valid_input = true;
		while (!fin.eof() && valid_input) {
			vector<GraphClass> g_list(BUFFER_SIZE);
			unsigned i = 0;
			while (i < BUFFER_SIZE && !fin.eof() && valid_input) {
				mpData->SetGraphFromFile(fin, g_list[i]);
				if (g_list[i].IsEmpty())
					valid_input = false;
				else
					i++;
			}
#pragma omp parallel for schedule(dynamic,100) ordered
			for (unsigned j = 0; j < i; j++) {
#pragma omp ordered
				{
					SVector x;
					mpData->mKernel.GenerateFeatureVector(g_list[j], x);
					if (mpParameters->mBinaryFormat)
						x.save(ofs);
					else
						ofs << x;
					pb.Count();
					instance_counter++;
				}
			}
		}
	}

};

//------------------------------------------------------------------------------------------------------------------------
class FeaturePartManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
public:
	FeaturePartManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mpData->mKernel.ParametersSetup();
	}

	void Exec() {
		InputOutputManager();
	}

	void InputOutputManager() {
		//output
		//compose output filename
		string output_filename = mpParameters->mInputDataFileName;
		output_filename += ".feature_part";
		output_filename += mpParameters->mSuffix;

		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.16: Cannot open file: " + output_filename);

		//input
		igzstream fin;
		fin.open(mpParameters->mInputDataFileName.c_str());
		if (!fin)
			throw range_error("ERROR2.11: Cannot open file: " + mpParameters->mInputDataFileName);
		//perform online action
		if (!mpParameters->mMinimalOutput)
			cout << "Processing I/O file: " << mpParameters->mInputDataFileName << endl;
		Main(fin, ofs);
		if (!mpParameters->mMinimalOutput)
			cout << "Result saved in file " + output_filename << endl;
	}

	void Main(istream& fin, ostream& ofs) {
		ProgressBar pb;
		unsigned instance_counter = 0;
		bool valid_input = true;
		while (!fin.eof() && valid_input) {
			vector<GraphClass> g_list(BUFFER_SIZE);
			unsigned i = 0;
			while (i < BUFFER_SIZE && !fin.eof() && valid_input) {
				mpData->SetGraphFromFile(fin, g_list[i]);
				if (g_list[i].IsEmpty())
					valid_input = false;
				else
					i++;
			}
#pragma omp parallel for schedule(dynamic,100) ordered
			for (unsigned j = 0; j < i; j++) {
#pragma omp ordered
				{
					vector<SVector> graph_vertex_vector_list;
					mpData->mKernel.GenerateVertexFeatureVector(g_list[j], graph_vertex_vector_list);
					unsigned size = mpParameters->mGraphType == "DIRECTED" ? graph_vertex_vector_list.size() / 2 : graph_vertex_vector_list.size();
					for (unsigned vertex_id = 0; vertex_id < size; ++vertex_id) {
						SVector x = graph_vertex_vector_list[vertex_id];
						if (mpParameters->mGraphType == "DIRECTED")
							x.add(graph_vertex_vector_list[vertex_id + size]);
						ofs << instance_counter << " " << vertex_id << " " << x;
						pb.Count();
						instance_counter++;
					}
				}
			}
		}
	}

};

//------------------------------------------------------------------------------------------------------------------------
class FeatureScaledManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	StochasticGradientDescentSupportVectorMachine mSGDSVM;

public:
	FeatureScaledManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mpData->mKernel.ParametersSetup();
		mSGDSVM.Init(mpParameters, mpData);
	}

	void LoadModel() {
		string filename = mpParameters->mModelFileName + mpParameters->mSuffix;
		ifstream ifs;
		ifs.open(filename.c_str());
		if (!ifs)
			throw range_error("ERROR2.23: Cannot open file:" + filename);
		if (!mpParameters->mMinimalOutput)
			cout << endl << "Loading model file: " << filename << endl;
		mSGDSVM.Load(ifs);
		if (!mpParameters->mMinimalOutput) {
			mSGDSVM.OutputModelInfo();
			cout << endl;
		}
	}

	void Exec() {
		LoadModel();
		InputOutputManager();
	}

	void InputOutputManager() {
		//output
		//compose output filename
		string output_filename = mpParameters->mInputDataFileName;
		output_filename += ".feature_scaled";
		output_filename += mpParameters->mSuffix;

		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.16: Cannot open file: " + output_filename);

		//input
		igzstream fin;
		fin.open(mpParameters->mInputDataFileName.c_str());
		if (!fin)
			throw range_error("ERROR2.11: Cannot open file: " + mpParameters->mInputDataFileName);
		//perform online action
		if (!mpParameters->mMinimalOutput)
			cout << "Processing I/O file: " << mpParameters->mInputDataFileName << endl;
		Main(fin, ofs);
		if (!mpParameters->mMinimalOutput)
			cout << "Result saved in file " + output_filename << endl;
	}

	void Main(istream& fin, ostream& ofs) {
		ProgressBar pb;
		unsigned instance_counter = 0;
		bool valid_input = true;
		while (!fin.eof() && valid_input) {
			switch (mpParameters->mFileTypeCode) {
			case GRAPH:
#ifdef USEOBABEL
			case MOLECULAR_GRAPH:
#endif
				case SEQUENCE: {
				vector<GraphClass> g_list(BUFFER_SIZE);
				unsigned i = 0;
				while (i < BUFFER_SIZE && !fin.eof() && valid_input) {
					mpData->SetGraphFromFile(fin, g_list[i]);
					if (g_list[i].IsEmpty())
						valid_input = false;
					else
						i++;
				}
#pragma omp parallel for schedule(dynamic,100) ordered
				for (unsigned j = 0; j < i; j++) {
#pragma omp ordered
					{
						SVector x;
						mpData->mKernel.GenerateFeatureVector(g_list[j], x);
						mSGDSVM.VectorElementwiseProductWithModel(x);
						if (mpParameters->mBinaryFormat)
							x.save(ofs);
						else
							ofs << x;
						pb.Count();
						instance_counter++;
					}
				}
			}
				break;
			case SPARSE_VECTOR: {
				SVector x;
				if (mpParameters->mBinaryFormat)
					mpData->SetVectorFromSparseVectorBinaryFile(fin, x);
				else
					mpData->SetVectorFromSparseVectorAsciiFile(fin, x);
				mSGDSVM.VectorElementwiseProductWithModel(x);
				if (mpParameters->mBinaryFormat)
					x.save(ofs);
				else
					ofs << x;
				pb.Count();
				instance_counter++;
			}
				break;
			default:
				throw range_error("ERROR2.45: file type not recognized: " + mpParameters->mFileType);
			}
		}
	}
}
;

//------------------------------------------------------------------------------------------------------------------------
class MinHashEncoder {
protected:
	Parameters* mpParameters;
	Data* mpData;
	vector<umap_uint_vec_uint> mInverseIndex;
	vector<vector<unsigned> > mMinHashCache;

public:
	MinHashEncoder() {
	}
	MinHashEncoder(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
	}

	void CacheReset() {
		mMinHashCache.clear();
		if (mpParameters->mNoMinHashCache == false)
			if (mpData->Size() > 0)
				mMinHashCache.resize(mpData->Size());

	}

	void ComputeInverseIndex() {
		CacheReset();
		cout << "Computing Inverse Index on " << mpData->mColIndexList.size() << " instances." << endl;
		cout << "Using " << mpParameters->mNumHashFunctions << " hash functions (with factor " << mpParameters->mNumRepeatsHashFunction << " for single minhash)" << endl;

		//init data structure
		mInverseIndex.clear();
		for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k)
			mInverseIndex.push_back(umap_uint_vec_uint());

		ProgressBar progress_bar;
#pragma omp parallel for schedule(dynamic)
		for (unsigned ii = 0; ii < mpData->mColIndexList.size(); ++ii) { //for every instance
			if (ii == 0) {
				int nthreads = omp_get_num_threads();
				cout << "Parallelization with " << nthreads << " threads" << endl;
			}

			unsigned i = mpData->mColIndexList[ii];
			vector<unsigned> signature = ComputeHashSignature(i); //compute the signature
#pragma omp critical
			{
				for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) { //for every hash value
					unsigned key = signature[k];
					if (key != MAXUNSIGNED && key != 0) { //if key is equal to markers for empty bins then skip insertion instance in data structure
						if (mInverseIndex[k].count(key) == 0) { //if this is the first time that an instance exhibits that specific value for that hash function, then store for the first time the reference to that instance
							vector<unsigned> tmp;
							tmp.push_back(i);
							mInverseIndex[k].insert(make_pair(key, tmp));
						} else { //Otherwise just append the instance to the list
							mInverseIndex[k][key].push_back(i);
						}
					}
				}
			}
			progress_bar.Count();
		}
	}

	inline vector<unsigned> ComputeHashSignature(unsigned aID) {
		if (mpParameters->mNoMinHashCache == false) {
			if (mMinHashCache[aID].size() > 0)
				return mMinHashCache[aID];
			else {
				vector<unsigned> signature = ComputeHashSignature(mpData->mVectorList[aID]);
				mMinHashCache[aID] = signature;
				return signature;
			}
		} else
			return ComputeHashSignature(mpData->mVectorList[aID]);
	}

	inline vector<unsigned> ComputeHashSignature(SVector& aX) {
		unsigned sub_hash_range = mpParameters->mNumHashFunctions / mpParameters->mNumRepeatsHashFunction;

		vector<unsigned> signature(mpParameters->mNumHashFunctions);
		//init with MAXUNSIGNED
		for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k)
			signature[k] = MAXUNSIGNED;

		//prepare a vector containing the signature as the k min values
		//for each element of the sparse vector
		for (unsigned f = 0; f < (unsigned) aX.sparse_size(); ++f) {
			unsigned feature_id = aX.extract_component(f).first;
			//for each sub_hash
			for (unsigned l = 1; l <= mpParameters->mNumRepeatsHashFunction; ++l) {
				unsigned key = IntHash(feature_id, MAXUNSIGNED, l);
				for (unsigned kk = 0; kk < sub_hash_range; ++kk) { //for all k values
					unsigned lower_bound = MAXUNSIGNED / sub_hash_range * kk;
					unsigned upper_bound = MAXUNSIGNED / sub_hash_range * (kk + 1);

					if (key >= lower_bound && key < upper_bound) { //if we are in the k-th slot
						unsigned signature_feature = kk + (l - 1) * sub_hash_range;
						if (key < signature[signature_feature]) //keep the min hash within the slot
							signature[signature_feature] = key;
					}
				}
			}
		}

		return signature;
	}

	vector<unsigned> ComputeApproximateNeighborhood(const vector<unsigned>& aSignature) {
		umap_uint_int neighborhood;
		for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) {
			unsigned hash_id = aSignature[k];
			if (hash_id != 0 && hash_id != MAXUNSIGNED) { //do not consider buckets corresponding to null bins
				unsigned collision_size = mInverseIndex[k][hash_id].size();

				if (collision_size < mpParameters->mMaxSizeBin * mpData->Size()) {
					//fill neighborhood set counting number of occurrences
					for (vector<unsigned>::iterator it = mInverseIndex[k][hash_id].begin(); it != mInverseIndex[k][hash_id].end(); ++it) {
						unsigned instance_id = *it;
						if (neighborhood.count(instance_id) > 0)
							neighborhood[instance_id]++;
						else
							neighborhood[instance_id] = 1;
					}
				} else {
					mInverseIndex[k].erase(hash_id);
				}
			}
		}
		return TrimNeighborhood(neighborhood);
	}

	vector<unsigned> TrimNeighborhood(umap_uint_int& aNeighborhood) {
		const int MIN_BINS_IN_COMMON = 2; //Minimum number of bins that two instances have to have in common in order to be considered similar
		//given a list of neighbors with an associated occurrences count, return only a fraction of the highest count ones
		vector<unsigned> neighborhood_list;
		if (mpParameters->mEccessNeighbourSizeFactor > 0) {
			//sort by num occurrences
			vector<pair<int, unsigned> > count_list;
			for (umap_uint_int::const_iterator it = aNeighborhood.begin(); it != aNeighborhood.end(); ++it) {
				unsigned id = it->first;
				int count = it->second;
				if (count >= MIN_BINS_IN_COMMON) //NOTE: consider instances that have at least MIN_BINS_IN_COMMON
					count_list.push_back(make_pair(-count, id)); //NOTE:-count to sort from highest to lowest
			}
			sort(count_list.begin(), count_list.end());
			unsigned effective_size = min((unsigned) count_list.size(), (unsigned) (mpParameters->mEccessNeighbourSizeFactor * mpParameters->mNumNearestNeighbors));
			for (unsigned i = 0; i < effective_size; ++i)
				neighborhood_list.push_back(count_list[i].second);
		} else { //if mEccessNeighbourSizeFactor==0 then just consider all the ids in the approximate neighborhood
			for (umap_uint_int::const_iterator it = aNeighborhood.begin(); it != aNeighborhood.end(); ++it) {
				neighborhood_list.push_back(it->first);
			}
		}
		return neighborhood_list;
	}
}
;

//------------------------------------------------------------------------------------------------------------------------
class MinHashManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	MinHashEncoder mMinHashEncoder;

public:
	MinHashManager() {
	}
	MinHashManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mpData->mKernel.ParametersSetup();
		mMinHashEncoder.Init(apParameters, apData);
	}

	void Exec() {
		InputOutputManager();
	}

	void InputOutputManager() {
		//output
		//compose output filename
		string output_filename = mpParameters->mInputDataFileName;
		output_filename += ".min_hash_feature";
		output_filename += mpParameters->mSuffix;

		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.167: Cannot open file: " + output_filename);

		//input
		igzstream fin;
		fin.open(mpParameters->mInputDataFileName.c_str());
		if (!fin)
			throw range_error("ERROR2.117: Cannot open file: " + mpParameters->mInputDataFileName);
		//perform online action
		if (!mpParameters->mMinimalOutput)
			cout << "Processing I/O file: " << mpParameters->mInputDataFileName << endl;
		Main(fin, ofs);
		if (!mpParameters->mMinimalOutput)
			cout << "Result saved in file " + output_filename << endl;
	}

	void Main(istream& fin, ostream& ofs) {
		ProgressBar pb;
		bool valid_input = true;
		while (!fin.eof() && valid_input) {

			switch (mpParameters->mFileTypeCode) {
			case GRAPH:
#ifdef USEOBABEL
			case MOLECULAR_GRAPH:
#endif
				case SEQUENCE: {
				GraphClass g;
				mpData->SetGraphFromFile(fin, g);
				if (!g.IsEmpty()) {
					SVector x;
					mpData->mKernel.GenerateFeatureVector(g, x);
					vector<unsigned> signature = mMinHashEncoder.ComputeHashSignature(x); //compute the signature
					for (unsigned k = 0; k < signature.size(); k++) {
						ofs << signature[k] << " ";
					}
					ofs << endl;
					pb.Count();
				} else
					valid_input = false;
			}
				break;
			default:
				throw range_error("ERROR2.45: file type not recognized: " + mpParameters->mFileType);
			}
		}
	}
};

//------------------------------------------------------------------------------------------------------------------------
class NearestNeighborManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	MinHashEncoder mMinHashEncoder;

public:
	NearestNeighborManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mMinHashEncoder.Init(apParameters, apData);
	}

	void Load() {
		mpData->LoadData();
	}

	void Exec() {
		Load();
		OutputManager();
	}

	void OutputManager() {
		cout << SEP << endl << "K-Nearest-Neighbor phase" << endl << SEP << endl;
		ProgressBar pb;
		pb.Count();

		string output_filename_knn = mpParameters->mInputDataFileName + ".knn" + mpParameters->mSuffix;
		ofstream ofs_knn;
		ofs_knn.open(output_filename_knn.c_str());
		if (!ofs_knn)
			throw range_error("ERROR2.16: Cannot open file:" + output_filename_knn);

		string output_filename_k = mpParameters->mInputDataFileName + ".knn_kernel_value" + mpParameters->mSuffix;
		ofstream ofs_k;
		ofs_k.open(output_filename_k.c_str());
		if (!ofs_k)
			throw range_error("ERROR2.16: Cannot open file:" + output_filename_k);

		string output_filename_t = mpParameters->mInputDataFileName + ".knn_target_value" + mpParameters->mSuffix;
		ofstream ofs_t;
		if (mpParameters->mTargetFileName != "") {
			ofs_t.open(output_filename_t.c_str());
			if (!ofs_t)
				throw range_error("ERROR2.16: Cannot open file:" + output_filename_t);
		}

		Main(ofs_knn, ofs_k, ofs_t);
		cout << "Nearest neighbor list saved in file " << output_filename_knn << endl;
		cout << "Nearest neighbor kernel value list saved in file " << output_filename_k << endl;
		if (mpParameters->mTargetFileName != "")
			cout << "Nearest neighbor target value list saved in file " << output_filename_t << endl;
	}

	void Main(ostream& out_knn, ostream& out_k, ostream& out_t) {
		if (mpData->mRowIndexList.size() == 0) {
			cout << "No row index list specified. Assuming all " << mpData->Size() << " row indices as valid." << endl;
			for (unsigned i = 0; i < mpData->Size(); ++i)
				mpData->mRowIndexList.push_back(i);
		}
		if (mpData->mColIndexList.size() == 0) {
			cout << "No col index list specified. Assuming all " << mpData->Size() << " col indices as valid." << endl;
			for (unsigned i = 0; i < mpData->Size(); ++i)
				mpData->mColIndexList.push_back(i);
		}

		if (mpParameters->mForceApproximate == true)
			mMinHashEncoder.ComputeInverseIndex();

		cout << "Computing " << mpParameters->mNumNearestNeighbors << " nearest neighbors for " << mpData->mRowIndexList.size() * mpData->mColIndexList.size() << " pairs of instances." << endl;

		ProgressBar ppb;
		for (unsigned i = 0; i < mpData->mRowIndexList.size(); ++i) {
			unsigned ii = mpData->mRowIndexList[i];
			vector<pair<double, unsigned> > neighbor_list;

			if (mpParameters->mForceApproximate == false) { //do full comparison of each row instance vs. each column instance and sort to find nearest
				vector<pair<double, unsigned> > init_neighbor_list(mpData->mColIndexList.size());
				neighbor_list = init_neighbor_list;
#pragma omp parallel for schedule(dynamic)
				for (unsigned j = 0; j < mpData->mColIndexList.size(); ++j) {
					if (j == 0 && i == 0) {
						int nthreads = omp_get_num_threads();
						cout << "Parallelization with " << nthreads << " threads" << endl;
					}
					unsigned jj = mpData->mColIndexList[j];
					double k_ii_jj = mpData->ComputeKernel(ii, jj);
					neighbor_list[j] = make_pair(-k_ii_jj, jj); //NOTE: use -k to sort in descending order of similarity
				}
				sort(neighbor_list.begin(), neighbor_list.end());
			} else { // extract approximate neighbors
				SVector& x = mpData->mVectorList[ii];
				vector<unsigned> signature = mMinHashEncoder.ComputeHashSignature(x);
				vector<unsigned> approximate_neighborhood = mMinHashEncoder.ComputeApproximateNeighborhood(signature);
				vector<pair<double, unsigned> > rank_list(approximate_neighborhood.size());
#pragma omp parallel for schedule(dynamic,10)
				for (unsigned j = 0; j < approximate_neighborhood.size(); ++j) {
					if (j == 0 && i == 0) {
						int nthreads = omp_get_num_threads();
						cout << "Parallelization with " << nthreads << " threads" << endl;
					}
					unsigned jj = approximate_neighborhood[j];
					double k = mpData->mKernel.ComputeKernel(mpData->mVectorList[jj], x);
					rank_list[j] = make_pair(-k, jj);
				}
				unsigned effective_size = min((unsigned) rank_list.size(), mpParameters->mNumNearestNeighbors);
				partial_sort(rank_list.begin(), rank_list.begin() + effective_size, rank_list.end());
				vector<pair<double, unsigned> > init_neighbor_list(effective_size);
				neighbor_list = init_neighbor_list;
				for (unsigned j = 0; j < effective_size; j++) {
					neighbor_list[j] = rank_list[j];
				}
			}
			ppb.Count();

			unsigned effective_size = min((unsigned) neighbor_list.size(), mpParameters->mNumNearestNeighbors);
			out_knn << ii << "    ";
			out_k << ii << "    ";
			if (mpParameters->mTargetFileName != "")
				out_t << ii << "    ";

			for (unsigned t = 0; t < effective_size; ++t) {
				double neighbor_kernel = -neighbor_list[t].first; //NOTE:revert to positive kernel
				unsigned neighbor_id = neighbor_list[t].second;
				out_knn << neighbor_id << " ";
				out_k << neighbor_kernel << " ";
				if (mpParameters->mTargetFileName != "")
					out_t << mpData->mTargetList[neighbor_id] << " ";
			}
			out_knn << endl;
			out_k << endl;
			if (mpParameters->mTargetFileName != "")
				out_t << endl;
		}
	}

};

//------------------------------------------------------------------------------------------------------------------------
class ClusterManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	MinHashEncoder mMinHashEncoder;
	vector<vector<unsigned> > mNeighbourhoodCache;
public:
	ClusterManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mpData->LoadData();
		mMinHashEncoder.Init(apParameters, apData);
		mMinHashEncoder.ComputeInverseIndex();

		if (mpParameters->mNoNeighborhoodCache == false)
			CacheReset();
	}

	void CacheReset() {
		mNeighbourhoodCache.clear();
		mNeighbourhoodCache.resize(mpData->Size());
	}

	void Exec() {
		ProgressBar pb;

		if (mpParameters->mClusterType == "DENSE_CENTERS") {
			string output_filename = mpParameters->mInputDataFileName;
			output_filename += ".fast_cluster" + mpParameters->mSuffix;
			ofstream ofs(output_filename.c_str());

			string output_filename_sim = mpParameters->mInputDataFileName;
			output_filename_sim += ".fast_cluster_inner_similarity" + mpParameters->mSuffix;
			ofstream ofs_sim(output_filename_sim.c_str());

			string output_filename_eccess = mpParameters->mInputDataFileName;
			output_filename_eccess += ".fast_cluster_eccess_neighbours" + mpParameters->mSuffix;
			ofstream ofs_ecc(output_filename_eccess.c_str());

			DenseCluster(ofs, ofs_sim, ofs_ecc);
			cout << endl << "Dense center cluster results written in file " << output_filename << endl;
			cout << endl << "Similarity matrix results written in file " << output_filename_sim << endl;
			cout << endl << "Eccess neighbours for dense centers written in file " << output_filename_eccess << endl;
		}

		if (mpParameters->mClusterType == "K_QUICK_SHIFT") {

			string output_filename = mpParameters->mInputDataFileName;
			output_filename += ".kquickshift_cluster";
			output_filename += mpParameters->mSuffix;
			ofstream ofs(output_filename.c_str());

			string output_filename_parent = mpParameters->mInputDataFileName;
			output_filename_parent += ".kquickshift_cluster_parent_relation";
			output_filename_parent += mpParameters->mSuffix;
			ofstream ofsp(output_filename_parent.c_str());

			string output_filename_ext = mpParameters->mInputDataFileName;
			output_filename_ext += ".kquickshift_cluster_neighbour_extension";
			output_filename_ext += mpParameters->mSuffix;
			ofstream ofsext(output_filename_ext.c_str());

			string output_filename_sim = mpParameters->mInputDataFileName;
			output_filename_sim += ".kquickshift_cluster_inner_similarity" + mpParameters->mSuffix;
			ofstream ofs_sim(output_filename_sim.c_str());

			KQuickShiftCluster(ofs, ofsp, ofsext, ofs_sim);
			cout << endl << "Clustering results written in file " << output_filename << endl;
			cout << endl << "Parent relations written in file " << output_filename_parent << endl;
			cout << endl << "Nearest Neighbor Extension Clustering results written in file " << output_filename_ext << endl;
			cout << endl << "Similarity results for each cluster written in file " << output_filename_sim << endl;
		}

		if (mpParameters->mClusterType == "APPROXIMATION_ACCURACY") {
			string output_filename = mpParameters->mInputDataFileName;
			output_filename += ".nn_accuracy";
			output_filename += mpParameters->mSuffix;
			ofstream ofs(output_filename.c_str());
			Accuracy(ofs);
			cout << endl << "Results written in file " << output_filename << endl;
		}
		pb.Count();
		cout << "Total clustering time:" << endl;
	}

	vector<unsigned> ComputeTrueSubNeighborhood(unsigned aID, vector<unsigned>& aApproximateNeighborhoodList) {
		//returns a subset of the approximate neighbours sorted by the true kernel function
		vector<unsigned> neighborhood_list;
		vector<pair<double, unsigned> > rank_list;
		for (unsigned i = 0; i < aApproximateNeighborhoodList.size(); ++i) {
			unsigned id_neighbour = aApproximateNeighborhoodList[i];
			double k = mpData->ComputeKernel(aID, id_neighbour);
			rank_list.push_back(make_pair(-k, id_neighbour));
		}
		unsigned effective_size = min((unsigned) rank_list.size(), mpParameters->mNumNearestNeighbors);
		partial_sort(rank_list.begin(), rank_list.begin() + effective_size, rank_list.end());

		for (unsigned j = 0; j < effective_size; j++) {
			neighborhood_list.push_back(rank_list[j].second);
		}
		return neighborhood_list;
	}

	vector<unsigned> ComputeTrueNeighborhood(unsigned aID) {
		const SVector& x = mpData->mVectorList[aID];
		return ComputeTrueNeighborhood(x);
	}

	vector<unsigned> ComputeTrueNeighborhood(const SVector& aX) {
		vector<pair<double, unsigned> > rank_list(mpData->Size());
#pragma omp parallel for schedule(dynamic)
		for (unsigned i = 0; i < mpData->Size(); ++i) {
			double k = mpData->mKernel.ComputeKernel(mpData->mVectorList[i], aX);
			rank_list[i] = make_pair(-k, i);
		}
		unsigned effective_size = min((unsigned) rank_list.size(), mpParameters->mNumNearestNeighbors);
		partial_sort(rank_list.begin(), rank_list.begin() + effective_size, rank_list.end());
		vector<unsigned> neighborhood_list;
		for (unsigned j = 0; j < effective_size; j++)
			neighborhood_list.push_back(rank_list[j].second);
		return neighborhood_list;
	}

	vector<unsigned> ComputeApproximateNeighborhood(unsigned aID) {
		vector<unsigned> approximate_true_neighborhood;
		if (mpParameters->mNoNeighborhoodCache == true) {
			vector<unsigned> signature = mMinHashEncoder.ComputeHashSignature(aID);
			vector<unsigned> approximate_neighborhood = mMinHashEncoder.ComputeApproximateNeighborhood(signature);
			approximate_true_neighborhood = ComputeTrueSubNeighborhood(aID, approximate_neighborhood);
		} else {
			if (mNeighbourhoodCache[aID].size() != 0) {
				approximate_true_neighborhood = mNeighbourhoodCache[aID];
			} else {
				vector<unsigned> signature = mMinHashEncoder.ComputeHashSignature(aID);
				vector<unsigned> approximate_neighborhood = mMinHashEncoder.ComputeApproximateNeighborhood(signature);
				approximate_true_neighborhood = ComputeTrueSubNeighborhood(aID, approximate_neighborhood);
				mNeighbourhoodCache[aID] = approximate_true_neighborhood;
			}
		}
		return approximate_true_neighborhood;
	}

	double ComputeApproximateDensity(unsigned aID) {
		double density = 0;
		vector<unsigned> approximate_neighborhood = ComputeApproximateNeighborhood(aID);
		density = ComputeDensity(aID, approximate_neighborhood);
		return density;
	}

	double ComputeDensity(unsigned aID, vector<unsigned>& aNeighborhood) {
		double density = 0;
		//compute kernel pairs between i and all elements in aApproximateNeighborhood
		for (unsigned j = 0; j < aNeighborhood.size(); j++) {
			unsigned u = aID;
			unsigned v = aNeighborhood[j];
			if (u != v) {
				double k_uv = mpData->ComputeKernel(u, v);
				if (mpParameters->mSharedNeighborhood) //if we use the shared neighborhood weighting than we multiply the similarity by a corrective factor given by the fraction of shared neighbours between the two instances
					k_uv *= ComputeSharedNeighborhoodSimilarity(u, v);
				density += k_uv;
			}
		}
		density = density / (aNeighborhood.size() - 1); //-1 as the similarity of an instance to itself (which is the closest of the elements in the neighborhood)  is excluded
		return density;
	}

	double ComputeDensityStandardDeviation(unsigned aID, double aAverageDensity, vector<unsigned>& aNeighborhood) {
		double variance = 0;
		//compute kernel pairs between i and all elements in aApproximateNeighborhood
		for (unsigned j = 0; j < aNeighborhood.size(); j++) {
			unsigned u = aID;
			unsigned v = aNeighborhood[j];
			if (u != v) {
				double k_uv = mpData->ComputeKernel(u, v);
				if (mpParameters->mSharedNeighborhood) //if we use the shared neighborhood weighting than we multiply the similarity by a corrective factor given by the fraction of shared neighbours between the two instances
					k_uv *= ComputeSharedNeighborhoodSimilarity(u, v);
				variance += (k_uv - aAverageDensity) * (k_uv - aAverageDensity);
			}
		}
		variance = variance / (aNeighborhood.size() - 1); //-1 as the similarity of an instance to itself (which is the closest of the elements in the neighborhood)  is excluded
		return sqrt(variance);
	}

	vector<unsigned> ComputeMinimallyOverlappingHighDensityCenterList(unsigned aSampleSize, double aFractionCenterScan, unsigned aMaxIntersectionSize) {
		//compute density estimate for random fraction of instances in dataset
		unsigned data_size = mpData->mRowIndexList.size() > 0 ? mpData->mRowIndexList.size() : mpData->Size();
		unsigned effective_size = floor(data_size * aFractionCenterScan);
		if (aFractionCenterScan == 1)
			cout << "Selecting all the " << mpData->mRowIndexList.size() << " instances" << endl;
		else
			cout << "Selecting a random " << aFractionCenterScan * 100 << " % of " << mpData->mRowIndexList.size() << " instances" << endl;
		//random selection of instances
		vector<unsigned> selected_index_list;
		//select either all the available centers
		if (aFractionCenterScan == 1) {
			selected_index_list.insert(selected_index_list.begin(), mpData->mRowIndexList.begin(), mpData->mRowIndexList.end());
		} else {
			//or select a random subset of candidate centers of size effective_size
			vector<unsigned> index_list = mpData->mRowIndexList;
			for (unsigned i = 0; i < data_size; ++i) {
				unsigned j = randomUnsigned(data_size);
				swap(index_list[i], index_list[j]);
			}
			for (unsigned i = 0; i < effective_size; ++i)
				selected_index_list.push_back(index_list[i]);
		}

		//compute density estimate
		vector<pair<double, unsigned> > density_list;
		ComputeApproximateDensity(selected_index_list, density_list);
		//select non overlapping centers in decreasing order of density
		sort(density_list.begin(), density_list.end());

		vector<unsigned> result;
		set<unsigned> active_neighborhood;
		{
			ProgressBar progress_bar(1);
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
		return result;
	}

	void ComputeApproximateDensity(vector<unsigned>& aSelectedIndexList, vector<pair<double, unsigned> >& oDensityList) {
		cout << "Computing approximate density information for random sample of " << aSelectedIndexList.size() << " instances" << endl; ////

		ProgressBar progress_bar;

		oDensityList.clear();
		oDensityList.resize(aSelectedIndexList.size());

#pragma omp parallel for schedule(dynamic)
		for (unsigned j = 0; j < aSelectedIndexList.size(); ++j) {
			if (j == 0) {
				int nthreads = omp_get_num_threads();
				cout << "Parallelization with " << nthreads << " threads" << endl;
			}
			unsigned i = aSelectedIndexList[j];
			double density = ComputeApproximateDensity(i);
			oDensityList[j] = make_pair(-density, i);
			progress_bar.Count();
		}
	}

	void DenseCluster(ostream& out, ostream& out_sim, ostream& out_ecc) {
		vector<unsigned> density_center_list;
		density_center_list = ComputeMinimallyOverlappingHighDensityCenterList(mpParameters->mSampleSize, mpParameters->mFractionCenterScan, mpParameters->mMaxIntersectionSize);
		cout << "Compute (approximate) neighborhood for selected " << density_center_list.size() << " cluster centers." << endl;
		{
			ProgressBar progress_bar(1);
			for (unsigned idc = 0; idc < density_center_list.size(); ++idc) {
				unsigned id = density_center_list[idc];
				vector<unsigned> neighborhood;
				if (mpParameters->mForceApproximate)
					neighborhood = ComputeApproximateNeighborhood(id);
				else
					neighborhood = ComputeTrueNeighborhood(id);
				progress_bar.Count();
				OutputDenseClusterNeighbors(out, neighborhood);
				OutputDenseClusterNeighborsSimilarity(out_sim, neighborhood);

				//output excess neighborhood
				if (mpParameters->mNoNeighborhoodCache == false)
					CacheReset();
				unsigned num_nearest_neighbours = mpParameters->mNumNearestNeighbors;
				mpParameters->mNumNearestNeighbors = mpParameters->mEccessNeighbourSizeFactor * mpParameters->mNumNearestNeighbors;
				vector<unsigned> neighborhood_eccess;
				neighborhood_eccess = ComputeApproximateNeighborhood(id);
				OutputDenseClusterNeighbors(out_ecc, neighborhood_eccess);
				mpParameters->mNumNearestNeighbors = num_nearest_neighbours;
			}
		}
	}

	void OutputDenseClusterNeighbors(ostream& out, vector<unsigned>& aNeighborhood) {
		for (unsigned i = 0; i < aNeighborhood.size(); i++) {
			unsigned nid = aNeighborhood[i];
			out << nid << " ";
		}
		out << endl;
	}

	void OutputDenseClusterNeighborsSimilarity(ostream& out, vector<unsigned>& aNeighborhood) {
		for (unsigned i = 0; i < aNeighborhood.size(); i++) {
			unsigned nid = aNeighborhood[i];
			for (unsigned j = i + 1; j < aNeighborhood.size(); j++) {
				unsigned njd = aNeighborhood[j];
				double k = mpData->ComputeKernel(nid, njd);
				out << nid << ":" << njd << ":" << k << " ";
			}
			out << " ";
		}
		out << endl;
	}

	pair<double, double> ComputeLocalCentrality(unsigned aID, vector<vector<unsigned> >& aApproximateNeighbourhoodList) {
		double centrality = ComputeDensity(aID, aApproximateNeighbourhoodList[aID]);
		double std = ComputeDensityStandardDeviation(aID, centrality, aApproximateNeighbourhoodList[aID]);
		return make_pair(centrality, std);
	}

	void KQuickShiftCluster(ostream& out, ostream& out_parent, ostream& out_ext, ostream& out_sim) {
		//compute approximate neighborhood and cache it
		vector<vector<unsigned> > approximate_neighborhood_list(mpData->Size());

		cout << "Compute approximate neighborhood for " << mpData->mRowIndexList.size() << " instances." << endl;
		{
			ProgressBar progress_bar;
#pragma omp parallel for schedule(dynamic) ordered
			for (unsigned ii = 0; ii < mpData->mRowIndexList.size(); ii++) {
				if (ii == 0) {
					int nthreads = omp_get_num_threads();
					cout << "Parallelization with " << nthreads << " threads" << endl;
				}
				unsigned i = mpData->mRowIndexList[ii];
				vector<unsigned> approximate_neighborhood = ComputeApproximateNeighborhood(i);
#pragma omp ordered
				{
					approximate_neighborhood_list[i] = approximate_neighborhood;
					progress_bar.Count();
				}
			}
		}

		//compute local centrality and cache it
		vector<double> local_centrality_list(mpData->Size());
		vector<double> local_centrality_std_list(mpData->Size());
		cout << "Compute local centality for " << mpData->mRowIndexList.size() << " instances." << endl;
		{
			ProgressBar progress_bar;
#pragma omp parallel for schedule(dynamic) ordered
			for (unsigned ii = 0; ii < mpData->mRowIndexList.size(); ii++) {
				unsigned i = mpData->mRowIndexList[ii];
				pair<double, double> stats = ComputeLocalCentrality(i, approximate_neighborhood_list);
#pragma omp ordered
				{
					double local_centralty = stats.first;
					double local_centrality_std = stats.second;
					local_centrality_list[i] = local_centralty;
					local_centrality_std_list[i] = local_centrality_std;
					progress_bar.Count();
				}
			}
		}

		//find closest more central neighbor and store parent pointer
		vector<int> parent_pointer_list(mpData->Size(), -1);
		cout << "Compute parent pointer for " << mpData->mRowIndexList.size() << " instances." << endl;
		{
			ProgressBar progress_bar;
			for (unsigned ii = 0; ii < mpData->mRowIndexList.size(); ii++) {
				unsigned i = mpData->mRowIndexList[ii];
				unsigned parent_id = i;
				double local_centrality = local_centrality_list[i];
				for (unsigned j = 0; j < approximate_neighborhood_list[i].size(); j++) {
					unsigned neighbour_id = approximate_neighborhood_list[i][j];
					if (local_centrality_list[neighbour_id] > local_centrality) {
						parent_id = neighbour_id;
						break;
					}
				}
				parent_pointer_list[i] = parent_id;
				out_parent << i << " " << parent_id << endl;

				progress_bar.Count();
			}
		}

		//invert parent  relationship and compute child relationship
		vector<vector<unsigned> > child_pointer_list(mpData->Size());
		for (unsigned ii = 0; ii < mpData->mRowIndexList.size(); ii++) {
			unsigned i = mpData->mRowIndexList[ii];
			int parent_id = parent_pointer_list[i];
			if (parent_id != -1) //if it is not a blacklisted parent
				if (parent_id != (int) i) { //if it is not a pointer to self
					child_pointer_list[parent_id].push_back(i);
				}
		}

		//for every instance compare average neighbor similarity and remove all parent edges such that parent has similarity below x normalized value
		vector<unsigned> root_pointer_list;
		vector<int> filtered_parent_pointer_list(mpData->Size());

		for (unsigned ii = 0; ii < mpData->mRowIndexList.size(); ii++) {
			unsigned i = mpData->mRowIndexList[ii];
			int parent_id = parent_pointer_list[i];
			if (parent_id == -1) {
			} //skip any blacklisted instance
			else {
				double average = local_centrality_list[i];
				double standard_deviation = local_centrality_std_list[i];
				double k_parent = mpData->ComputeKernel(i, parent_id);

				if (parent_id == (int) i) { //self loops (i.e. roots) are preserved
					filtered_parent_pointer_list[i] = i;
					root_pointer_list.push_back(i);
				} else if (child_pointer_list.size() == 0) { //leaf are preserved if the parent similarity is compatible with the average similarity around the parent
					double parent_average = local_centrality_list[parent_id];
					double parent_standard_deviation = local_centrality_std_list[parent_id];
					if (fabs(k_parent - parent_average) / parent_standard_deviation > 1 / mpParameters->mClusterThreshold) {
						filtered_parent_pointer_list[i] = parent_id;
					} else {
						filtered_parent_pointer_list[i] = i;
						root_pointer_list.push_back(i);
					}
				} else if (fabs(k_parent - average) / standard_deviation > 1 / mpParameters->mClusterThreshold) {
					filtered_parent_pointer_list[i] = parent_id;
				} else { //else make a new root
					filtered_parent_pointer_list[i] = i;
					root_pointer_list.push_back(i);
				}
			}
		}

//invert parent  relationship and compute filtered child relationship
		vector<vector<unsigned> > filtered_child_pointer_list(mpData->Size());
		for (unsigned i = 0; i < mpData->Size(); i++) {
			int parent_id = filtered_parent_pointer_list[i];
			if (parent_id != -1) //if it is not a blacklisted parent
				if (parent_id != (int) i) //if it is not a pointer to self
					filtered_child_pointer_list[parent_id].push_back(i);
		}

//for each root visit the dominated tree and output all instances as members of the cluster
		for (unsigned i = 0; i < root_pointer_list.size(); ++i) {
			unsigned root_id = root_pointer_list[i];
			vector<unsigned> dominated_tree;
			TreeVisit(root_id, filtered_child_pointer_list, dominated_tree);
			for (unsigned j = 0; j < dominated_tree.size(); ++j) {
				out << dominated_tree[j] << " ";
			}
			out << endl;
		}

//for each elements in the kquickshift clusters, find knn and output them
		for (unsigned i = 0; i < root_pointer_list.size(); ++i) {
			set<unsigned> cluster_id_list;

			unsigned root_id = root_pointer_list[i];
			cluster_id_list.insert(approximate_neighborhood_list[root_id].begin(), approximate_neighborhood_list[root_id].end());

			vector<unsigned> dominated_tree;
			TreeVisit(root_id, filtered_child_pointer_list, dominated_tree);
			for (unsigned j = 0; j < dominated_tree.size(); ++j) {
				unsigned id = dominated_tree[j];
				cluster_id_list.insert(approximate_neighborhood_list[id].begin(), approximate_neighborhood_list[id].end());
			}
			for (set<unsigned>::iterator it = cluster_id_list.begin(); it != cluster_id_list.end(); ++it)
				out_ext << (*it) << " ";
			out_ext << endl;
		}

//for each elements in the kquickshift clusters, output their pairwise similarity
		for (unsigned i = 0; i < root_pointer_list.size(); ++i) {
			unsigned root_id = root_pointer_list[i];
			vector<unsigned> dominated_tree;
			TreeVisit(root_id, filtered_child_pointer_list, dominated_tree);
			if (dominated_tree.size() == 1) {
				out_sim << root_id << ":" << root_id << ":1" << endl;
			} else {
				for (unsigned j = 0; j < dominated_tree.size(); ++j) {
					unsigned idx = dominated_tree[j];
					for (unsigned z = j + 1; z < dominated_tree.size(); ++z) {
						unsigned idz = dominated_tree[z];
						double k = mpData->ComputeKernel(idx, idz);
						out_sim << idx << ":" << idz << ":" << k << " ";
					}
				}
				out_sim << endl;
			}
		}

//output also the average cluster similarity
	}

	void TreeVisit(unsigned aID, vector<vector<unsigned> >& aChildPointerList, vector<unsigned>& oDominatedTree) {
		oDominatedTree.push_back(aID);
		for (unsigned i = 0; i < aChildPointerList[aID].size(); ++i) {
			unsigned child_id = aChildPointerList[aID][i];
			TreeVisit(child_id, aChildPointerList, oDominatedTree);
		}
	}

	void Accuracy(ostream& out) {
		cout << "Compute neighbourhood accuracy" << endl; ////
		ProgressBar progress_bar;
		double cum = 0;
		unsigned effective_neighbourhood_size = min(mpParameters->mNumNearestNeighbors, (unsigned) mpData->Size());
		for (unsigned u = 0; u < mpData->Size(); ++u) {
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
		cout << endl << "Accuracy: " << cum / mpData->Size() << endl;
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
};

//------------------------------------------------------------------------------------------------------------------------
class SemiSupervisedManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
	MinHashEncoder mMinHashEncoder;

public:
	SemiSupervisedManager(Parameters* apParameters, Data* apData) {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData) {
		mpParameters = apParameters;
		mpData = apData;
		mMinHashEncoder.Init(apParameters, apData);
	}

	void Load() {
		mpData->LoadData();
	}

	void Exec() {
		Load();
		OutputManager();
	}

	void OutputManager() {
		cout << SEP << endl << "Semisupervised phase" << endl << SEP << endl;
		ProgressBar pb;
		pb.Count();

		string output_filename = mpParameters->mTargetFileName + ".semisupervised" + mpParameters->mSuffix;
		ofstream ofs;
		ofs.open(output_filename.c_str());
		if (!ofs)
			throw range_error("ERROR2.16: Cannot open file:" + output_filename);

		Main(ofs);
		cout << "Target value list saved in file " << output_filename << endl;
	}

	void Main(ostream& out) {
		if (mpParameters->mForceApproximate == true)
			mMinHashEncoder.ComputeInverseIndex();

		vector<SVector> neighbor_matrix;
		{
			cout << "Extracting nearest neighbor information" << endl;
			ProgressBar ppb;
			for (unsigned ii = 0; ii < mpData->Size(); ++ii) {
				vector<pair<double, unsigned> > neighbor_list;

				if (mpParameters->mForceApproximate == false) { //do full comparison of each row instance vs. each column instance and sort to find nearest
					vector<pair<double, unsigned> > init_neighbor_list(mpData->Size());
					neighbor_list = init_neighbor_list;
					for (unsigned jj = 0; jj < mpData->Size(); ++jj) {
						double k_ii_jj = mpData->ComputeKernel(ii, jj);
						neighbor_list[jj] = make_pair(-k_ii_jj, jj); //NOTE: use -k to sort in descending order of similarity
					}
					unsigned effective_size = min((unsigned) neighbor_list.size(), mpParameters->mNumNearestNeighbors);
					partial_sort(neighbor_list.begin(), neighbor_list.begin() + effective_size, neighbor_list.end());
				} else { // extract approximate neighbors
					vector<unsigned> signature = mMinHashEncoder.ComputeHashSignature(ii);
					vector<unsigned> approximate_neighborhood = mMinHashEncoder.ComputeApproximateNeighborhood(signature);
					vector<pair<double, unsigned> > rank_list;
					for (unsigned j = 0; j < approximate_neighborhood.size(); ++j) {
						unsigned jj = approximate_neighborhood[j];
						double k = mpData->ComputeKernel(jj, ii);
						rank_list.push_back(make_pair(-k, jj));
					}
					unsigned effective_size = min((unsigned) rank_list.size(), mpParameters->mNumNearestNeighbors);
					partial_sort(rank_list.begin(), rank_list.begin() + effective_size, rank_list.end());
					for (unsigned j = 0; j < effective_size; j++)
						neighbor_list.push_back(rank_list[j]);
				}
				ppb.Count();

				unsigned effective_size = min((unsigned) neighbor_list.size(), mpParameters->mNumNearestNeighbors);

				SVector effective_neighbor_list;
				double sum = 0;
				for (unsigned t = 0; t < effective_size; ++t) {
					double neighbor_kernel = -neighbor_list[t].first; //NOTE:revert to positive kernel
					unsigned neighbor_id = neighbor_list[t].second;
					if (neighbor_id != ii) { //avoid self loops
						effective_neighbor_list.set((int) neighbor_id, neighbor_kernel);
						sum += neighbor_kernel;
					}
				}
				//normalize
				effective_neighbor_list.scale(1 / sum);
				neighbor_matrix.push_back(effective_neighbor_list);
			}
		}

		//iterate activation spreading
		FVector y(mpData->mTargetList.size());
		//init y with target
		for (unsigned i = 0; i < mpData->mTargetList.size(); ++i)
			y.set(i, mpData->mTargetList[i]);
		//init f
		FVector f(mpData->mTargetList.size());
		f = y;
		//iterate
		{
			cout << "Spreading phase." << endl;
			ProgressBar ppb(1);
			for (unsigned itera = 0; itera < mpParameters->mSemiSupervisedNumIterations; ++itera) {
				ppb.Count();
				//spread
				FVector fprime(mpData->mTargetList.size());
				for (unsigned i = 0; i < neighbor_matrix.size(); ++i) {
					double val = dot(neighbor_matrix[i], f);
					fprime.set(i, val);
				}
				f = combine(fprime, mpParameters->mSemiSupervisedAlpha, y, (1 - mpParameters->mSemiSupervisedAlpha));
			}
		}
		//save
		for (int i = 0; i < f.size(); ++i) {
			out << f.get(i) << endl;
		}
	}

};

//------------------------------------------------------------------------------------------------------------------------
class Dispatcher {
protected:
	Parameters mParameters;
	Data mData;

public:
	Dispatcher() {
	}

	void Init(int argc, const char **argv) {
		mParameters.Init(argc, argv);
		srand(mParameters.mRandomSeed);
		mData.Init(&mParameters);
	}

	void Exec() {
		ProgressBar pb;

		if (!mParameters.mMinimalOutput)
			cout << SEP << endl << PROG_CREDIT << endl << SEP << endl;

		switch (mParameters.mActionCode) {
		case TRAIN: {
			StochasticGradientDescentSupportVectorMachineManager sgdsvm_manager(&mParameters, &mData);
			sgdsvm_manager.LoadTarget();
			sgdsvm_manager.LoadData();
			sgdsvm_manager.Train();
		}
			break;
		case CROSS_VALIDATION: {
			CrossValidationManager cross_validation_manager(&mParameters, &mData);
			cross_validation_manager.Exec();
		}
			break;
		case BIAS_VARIANCE_DECOMPOSITION: {
			BiasVarianceDecompositionManager bias_variance_decomposition_manager(&mParameters, &mData);
			bias_variance_decomposition_manager.Exec();
		}
			break;
		case PARAMETERS_OPTIMIZATION: {
			ParametersOptimizationManager parameters_optimization_manager(&mParameters, &mData);
			parameters_optimization_manager.Exec();
		}
			break;
		case LEARNING_CURVE: {
			LearningCurveManager learning_curve_manager(&mParameters, &mData);
			learning_curve_manager.Exec();
		}
			break;
		case TARGET_ALIGNMENT: {
			TargetAlignmentManager target_alignment_manager(&mParameters, &mData);
			target_alignment_manager.Exec();
		}
			break;
		case MATRIX: {
			GramMatrixManager gram_matrix_manager(&mParameters, &mData);
			gram_matrix_manager.Exec();
		}
			break;
		case NEAREST_NEIGHBOR: {
			NearestNeighborManager nearest_neighbor_manager(&mParameters, &mData);
			nearest_neighbor_manager.Exec();
		}
			break;
		case EMBED: {
			EmbedManager embed_manager(&mParameters, &mData);
			embed_manager.Exec();
		}
			break;
		case TEST: {
			TestManager test_manager(&mParameters, &mData);
			test_manager.Exec();
		}
			break;
		case TEST_PART: {
			TestPartManager test_part_manager(&mParameters, &mData);
			test_part_manager.Exec();
		}
			break;
		case FEATURE: {
			FeatureManager feature_manager(&mParameters, &mData);
			feature_manager.Exec();
		}
			break;
		case FEATURE_PART: {
			FeaturePartManager feature_part_manager(&mParameters, &mData);
			feature_part_manager.Exec();
		}
			break;
		case FEATURE_SCALED: {
			FeatureScaledManager feature_scaled_manager(&mParameters, &mData);
			feature_scaled_manager.Exec();
		}
			break;
		case CLUSTER: {
			ClusterManager cluster_manager(&mParameters, &mData);
			cluster_manager.Exec();
		}
			break;
		case MIN_HASH: {
			MinHashManager min_hash_manager(&mParameters, &mData);
			min_hash_manager.Exec();
		}
			break;
		case SEMI_SUPERVISED: {
			SemiSupervisedManager semi_supervised_manager(&mParameters, &mData);
			semi_supervised_manager.Exec();
		}
			break;
		default:
			throw range_error("ERROR2.2: Unknown action parameter: " + mParameters.mAction);
		}
		pb.Count();
		cout << "Total run-time:" << endl;
	}
};

//------------------------------------------------------------------------------------------------------------------------
int main(int argc, const char **argv) {
	try {
		Dispatcher dispatcher;
		dispatcher.Init(argc, argv);
		dispatcher.Exec();
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
	return 0;
}

//TODO: make unique output file manager that uses prefix and knows name of class
//TODO: make a better data manager that can be used also for in-line processing with an external callable function inbetween
//TODO: make parameter optimization
//TODO: make a principled object oriented organization of the code, so that the learning algorithm is a module, the instance representation is abstracted (i.e. accessible only through properties)
