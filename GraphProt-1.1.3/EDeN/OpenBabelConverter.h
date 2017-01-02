#ifndef CONVERTS2G_H_
#define CONVERTS2G_H_

#include "Utility.h"
#include "BaseGraphClass.h"
#include "GraphClass.h"

#ifdef USEOBABEL
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#endif

using namespace std;


#ifdef USEOBABEL
class OpenBabelConverter
{
  OpenBabel::OBConversion obconversion;
  bool	mInitStatus;
public:
  void ConvertOpenBabelFormatToGraph(istream* pIn, GraphClass& oG, string aFormat);
  void ConvertOpenBabelFormatToGraph(string aFileName, GraphClass& oG, string aFormat);
  void Convert(OpenBabel::OBMol& aObMol, GraphClass& oG);
  ~OpenBabelConverter() {};
};
#endif
#ifndef USEOBABEL
class OpenBabelConverter
{
public:
  void ConvertOpenBabelFormatToGraph(istream* pIn, GraphClass& oG, string aFormat){}
  void ConvertOpenBabelFormatToGraph(string aFileName, GraphClass& oG, string aFormat){}
};
#endif




#endif
