#include "Utility.h"

//------------------------------------------------------------------------------------------------------------------------
ParameterType::ParameterType() :
		mShortSwitch(""), mLongSwitch(""), mShortDescription(""), mLongDescription(""), mTypeCode(STRING), mValue(""), mMinValue(""), mMaxValue(""), mNumStepsValue(""), mExponentialStepIncrease(false), mIsSet(false) {
}

void ParameterType::Parse(vector<string>& aParameterList) {
	string shortopt = "-" + mShortSwitch;
	string longopt = "--" + mLongSwitch;
	string shortopt_min = "-" + mShortSwitch + "-min";
	string shortopt_max = "-" + mShortSwitch + "-max";
	string longopt_min = "-" + mLongSwitch + "-min";
	string longopt_max = "-" + mLongSwitch + "-max";
	string shortopt_numsteps = "-" + mShortSwitch + "-num_steps";
	string longopt_numsteps = "-" + mLongSwitch + "-num_steps";

	string str_value="";

	for (unsigned i = 0; i < aParameterList.size(); ++i) {
		if (aParameterList[i] == shortopt || aParameterList[i] == longopt || aParameterList[i] == shortopt_min || aParameterList[i] == shortopt_max || aParameterList[i] == longopt_min || aParameterList[i] == longopt_max || aParameterList[i] == shortopt_numsteps || aParameterList[i] == longopt_numsteps ) {
			mIsSet = true;
			if (mTypeCode != FLAG) {
				if (i + 1 >= aParameterList.size())
					throw range_error("ERROR: expected a value for option " + shortopt + "(" + longopt + ")");
				str_value = aParameterList[i + 1];
			}

			switch (mTypeCode) {
			case FLAG:
				str_value = "1";
				break;
			case REAL: {
				double value = strtod(str_value.c_str(), NULL);
				if (!value && str_value != "0")
					throw range_error("ERROR: Value for option " + shortopt + " (" + longopt + ") must be of type real");
			}
				break;
			case INTEGER: {
				int value = strtol(str_value.c_str(), NULL, 10);
				if (!value && str_value != "0")
					throw range_error("ERROR: Value for option " + shortopt + " (" + longopt + ") must be of type integer");
			}
				break;
			case POSITIVE_INTEGER: {
				unsigned value = strtoul(str_value.c_str(), NULL, 10);
				if (!value && str_value != "0")
					throw range_error("ERROR: Value for option " + shortopt + " (" + longopt + ") must be of type positive integer");
			}
				break;
			case STRING: {
			}
				break;
			case LIST: {
				bool test = false;
				for (unsigned i = 0; i < mCloseValuesList.size(); ++i) {
					if (mCloseValuesList[i] == str_value) {
						test = true;
						break;
					}
				}
				if (test == false) {
					string value_list = "";
					for (unsigned i = 0; i < mCloseValuesList.size() - 1; ++i)
						value_list += mCloseValuesList[i] + " | ";
					value_list += mCloseValuesList.back();
					throw range_error("ERROR: Value for option " + shortopt + " (" + longopt + ") must be one of: " + value_list + " instead of: " + mValue);
				}
			}
				break;
			default:
				throw range_error("ERROR: Option " + shortopt + " (" + longopt + ") must be set to an appropriate value");
			}
		}
		if (aParameterList[i] == shortopt || aParameterList[i] == longopt ){
			mValue=str_value;
		}
		if (aParameterList[i] == shortopt_min || aParameterList[i] == longopt_min ) {
			mMinValue=str_value;
		}
		if (aParameterList[i] == shortopt_max || aParameterList[i] == longopt_max ) {
			mMaxValue=str_value;
		}
		if (aParameterList[i] == shortopt_numsteps || aParameterList[i] == longopt_numsteps ) {
					mNumStepsValue=str_value;
		}
	}
}

void ParameterType::OutputCompact(ostream& out) const {
	if (mLongSwitch == "" && mShortSwitch == "")
		throw range_error("ERROR: Option does not have short nor long switch");
	if (mShortSwitch != "")
		out << "-" << mShortSwitch << " ";
	if (mLongSwitch != "")
		out << "--" << mLongSwitch << " ";
	if (mTypeCode == LIST) {
		string value_list = "";
		for (unsigned i = 0; i < mCloseValuesList.size() - 1; ++i)
			value_list += mCloseValuesList[i] + " | ";
		value_list += mCloseValuesList.back();
		out << value_list << " ";
	}
	if (mTypeCode != FLAG)
		if (mValue != "")
			out << "[" << mValue << "] ";
	if (mShortDescription != "")
		out << TAB << mShortDescription << " ";
	out << endl;
}

void ParameterType::OutputExtended(ostream& out) const {
	if (mLongSwitch == "" && mShortSwitch == "")
		throw range_error("ERROR: Option does not have short nor long switch");
	if (mShortSwitch != "")
		out << "-" << mShortSwitch << endl;
	if (mLongSwitch != "")
		out << "--" << mLongSwitch << endl;
	out << TAB << "Data type: ";
	switch (mTypeCode) {
	case FLAG:
		out << "flag" << endl;
		break;
	case REAL:
		out << "real number" << endl;
		break;
	case INTEGER:
		out << "integer number" << endl;
		break;
	case POSITIVE_INTEGER:
		out << "positive integer number" << endl;
		break;
	case STRING:
		out << "string" << endl;
		break;
	case LIST: {
		string value_list = "";
		for (unsigned i = 0; i < mCloseValuesList.size() - 1; ++i)
			value_list += mCloseValuesList[i] + " | ";
		value_list += mCloseValuesList.back();
		out << "One of the following " << mCloseValuesList.size() << " alternatives: " << value_list << endl;
	}
		break;
	default:
		throw range_error("ERROR4.1: Unrecognized parameter type code [" + mTypeCode);
	}
	if (mTypeCode != FLAG) {
		if (mValue != "")
			out << TAB << "Default: " << mValue << endl;
	}
	if (mLongDescription != "")
		out << TAB << "Description: " << mLongDescription << endl;
	else if (mShortDescription != "")
		out << TAB << "Description: " << mShortDescription << endl;
	out << endl;
}

//----------------------------------------------------------------------------------------------------------------------------
//hash functions

unsigned RSHash(const string& aString) {
	unsigned int b = 378551;
	unsigned int a = 63689;
	unsigned int hash = 0;
	for (std::size_t i = 0; i < aString.length(); i++) {
		hash = hash * a + aString[i];
		a = a * b;
	}
	return hash;
}

unsigned RSHash(const vector<unsigned>& aV) {
	unsigned int b = 378551;
	unsigned int a = 63689;
	unsigned int hash = 0;
	for (std::size_t i = 0; i < aV.size(); i++) {
		hash = hash * a + aV[i];
		a = a * b;
	}
	return hash;
}

unsigned APHash(const string& aString) {
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aString.length(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aString[i] * (hash >> 3)) : (~(((hash << 11) + aString[i]) ^ (hash >> 5)));
	}
	return hash;
}

unsigned APHash(const vector<unsigned>& aV) {
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aV.size(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aV[i] * (hash >> 3)) : (~(((hash << 11) + aV[i]) ^ (hash >> 5)));
	}
	return hash;
}

unsigned HashFunc(const string& aString, unsigned aBitMask) { //NOTE: extract the least significant bits from the hash
	return APHash(aString) & aBitMask;
}

unsigned HashFunc(const vector<unsigned>& aList, unsigned aBitMask) {
	return APHash(aList) & aBitMask;
}

//------------------------------------------------------------------------------------------------------------------------
TimerClass::TimerClass() :
		start_sec(time(NULL)), start(std::clock()) {
}
void TimerClass::Output() {
	//std::clock_t end=std::clock();
	//std::clock_t elapsed=end-start;
	std::time_t end_sec = time(NULL);
	double elapsed_sec = end_sec - start_sec;
	double elapsed_min = elapsed_sec / 60;
	double elapsed_hour = elapsed_min / 60;
	//CLOCKS_PER_SEC
	cout << "Elapsed time: h: " << floor(elapsed_hour) << " m: " << floor(elapsed_min) << " s: " << floor(elapsed_sec) << endl;
}

//------------------------------------------------------------------------------------------------------------------------
ProgressBar::ProgressBar(unsigned aStep) :
		mStep(aStep), mCounter(0) {
}
ProgressBar::~ProgressBar() {
	if (mCounter > 1) {
		cout << endl << "Counted " << mCounter << " times." << endl;
	}
	mTimer.Output();
}
void ProgressBar::Begin() {
	mCounter = 0;
}
void ProgressBar::Count() {
	mCounter++;
	if (mCounter % mStep == 0)
		cout << "." << flush;
	if (mCounter % 1000 == 0)
		cout << mCounter / (1000) << "K" << flush;
}
unsigned ProgressBar::End() {
	return mCounter;
}

//------------------------------------------------------------------------------------------------------------------------
ostream& operator<<(ostream& out, const VectorClass& aV) {
	aV.Output(out);
	return out;
}
VectorClass::VectorClass() {
}
VectorClass::VectorClass(unsigned aSize) {
	Init(aSize);
}
void VectorClass::operator=(const VectorClass& aVector) {
	mV = aVector.mV;
}
VectorClass::VectorClass(const VectorClass& aVector) {
	(*this) = aVector;
}
VectorClass::VectorClass(const vector<double>& aVector) {
	mV = aVector;
}
void VectorClass::Init(unsigned aSize) {
	mV.clear();
	for (unsigned i = 0; i < aSize; ++i)
		mV.push_back(-1);
}
void VectorClass::Import(const string& aFileName) {
	mV.clear();
	ifstream fin;
	fin.open(aFileName.c_str());
	if (!fin) {
		cerr << "Cannot open file: " << aFileName << endl;
		throw exception();
	}
	string line;

	//read size
	while (getline(fin, line)) {
		if (line != "") {
			stringstream ss;
			ss << line;
			while (ss.good()) {
				string value_str;
				ss >> value_str;
				if (value_str != "")
					mV.push_back(stream_cast<double>(value_str));
			}
		}
	}
	fin.close();
}
void VectorClass::Clear() {
	mV.clear();
}
unsigned VectorClass::Size() const {
	return mV.size();
}
ostream& VectorClass::Output(ostream& out) const {
	for (unsigned i = 0; i < Size(); i++)
		out << mV[i] << " ";
	return out;
}
void VectorClass::PushBack(double aValue) {
	mV.push_back(aValue);
}
double& VectorClass::operator[](unsigned i) {
	assert(i<Size());
	return mV[i];
}
double VectorClass::operator[](unsigned i) const {
	assert(i<Size());
	return mV[i];
}
double VectorClass::Sum() const {
	double avg = 0;
	for (unsigned i = 0; i < Size(); i++)
		avg += mV[i];
	return avg;
}
double VectorClass::Mean() const {
	return Sum() / Size();
}
double VectorClass::StandardDeviation() const {
	double avg = Mean();
	double sd = 0;
	for (unsigned i = 0; i < Size(); i++)
		sd += (mV[i] - avg) * (mV[i] - avg);
	sd = sd / (Size() - 1);
	sd = sqrt(sd);
	return sd;
}
double VectorClass::Order(double aOrder) const {
	vector<double> v(mV);
	sort(v.begin(), v.end());
	if (aOrder == 0)
		return v[0];
	else if (aOrder == 1)
		return v[v.size() - 1];
	else
		return v[(unsigned) ((double) Size() * aOrder)];
}
double VectorClass::Median() const {
	return Order(.5);
}
double VectorClass::MedianAbsoluteDifference() const {
	double median = Median();
	VectorClass v;
	for (unsigned i = 0; i < Size(); i++)
		v.PushBack(fabs(mV[i] - median));
	return v.Median();
}
double VectorClass::Min() const {
	return Order(0);
}
double VectorClass::Max() const {
	return Order(1);
}
VectorClass VectorClass::RemoveNulls() {
	VectorClass v;
	for (unsigned i = 0; i < Size(); i++)
		if (mV[i] != -1)
			v.PushBack(mV[i]);
	return v;
}
ostream& VectorClass::OutputStatistics(ostream& out) {
	VectorClass v = RemoveNulls();
	if (v.Size() > 0) {
		out << "num: " << v.Size() << " sum: " << v.Sum() << " avg: " << v.Mean() << " sd: " << v.StandardDeviation();
		out << " min: " << v.Min() << " Q.05: " << v.Order(.05) << " Q.25: " << v.Order(.25) << " Q.5: " << v.Median() << " Q.75: " << v.Order(.75) << " Q.95: " << v.Order(.95) << " max: " << v.Max();
	} else {
	}
	return out;
}

