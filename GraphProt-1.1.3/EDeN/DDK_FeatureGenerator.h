#ifndef CGRAPHPERCEPTRON_H
#define CGRAPHPERCEPTRON_H

/*
 * This program is free software; you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation; either version 3 of the License, or
  * (at your option) any later version.
  *
  * Written (W) 2013 Nicolò Navarin
  */

#include <stdio.h>
//#include <shogun/lib/common.h>
//#include <shogun/features/DotFeatures.h>
//#include <shogun/machine/LinearMachine.h>
//#include"CGraphFeatures.h"
//#include <boost/foreach.hpp>
#include "vectors.h"
//#include "bigDAG.h"
#include "NSPDK_FeatureGenerator.h"

using namespace std;


class DDkernel_FeatureGenerator :  public NSPDK_FeatureGenerator //public CLinearMachine,
 {
     public:
         DDkernel_FeatureGenerator();
         DDkernel_FeatureGenerator(const std::string& id);

       //  DDkernel_FeatureGenerator(CGraphFeatures* traindat, CLabels* trainlab);
         virtual ~DDkernel_FeatureGenerator();

       virtual  void generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList) ;
         void GenerateAbstractVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x);

         void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList);
       //  void svectorFeatures(CGraphFeatures* data,SVector& out_vector, int dag_h, int max_features);





         inline void set_lambda(double r)
         {
             mTreeLambda=r;
         }



         virtual void OutputParameters(ostream& out) const;

	unsigned getRadiusTwo() const {
		return mRadiusTwo;
	}

	void setRadiusTwo(unsigned radiusTwo) {
		mRadiusTwo = radiusTwo;
		cout<<"Set radius 2 "<<mRadiusTwo<<endl;

	}

     protected:
         double learn_rate;
         bool normalize;
         int32_t max_iter;
         double mTreeLambda;
         unsigned mRadiusTwo;
      //   double C;
 };


class DDkernel_FeatureGeneratorNew : public DDkernel_FeatureGenerator{
public:
    DDkernel_FeatureGeneratorNew(const std::string& id);
    void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList = vector<unsigned>());
    SVector&  generate_feature_vector_core(const GraphClass& aG , int dag_h);

};

class NSDDkernel_FeatureGenerator :  public DDkernel_FeatureGeneratorNew{
public:
	NSDDkernel_FeatureGenerator(const std::string& id);
	SVector&   generate_feature_vector_core(const GraphClass& aG , int dag_h);
  virtual  void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList = vector<unsigned>());

	SVector&  generate_vertex_feature_vector(unsigned i, const GraphClass& aG , int dag_h);
	virtual SVector&  generate_vertex_feature_vector_DIVIDED(unsigned i, const GraphClass& aG , int dag_h, vector<SVector>& features_grouped_by_depth);
	SVector&  generate_feature_vector_core_DISTANCE(const GraphClass& aG , int dag_h);

	SVector&  generate_vertex_feature_vector_DIVIDED_DISTANCE(unsigned i, const GraphClass& aG , int dag_h, vector<vector<SVector> > & features_grouped_by_depth);


};
class ANSDDkernel_FeatureGenerator :  public NSDDkernel_FeatureGenerator{
public:
	ANSDDkernel_FeatureGenerator(const std::string& id);
	void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList);

	//void generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList) ;
	void GenerateVertexFeatures(unsigned aRootVertexIndex, const GraphClass& aG, unsigned aRadius,unsigned aDistance, SVector& x);
	//void GenerateVertexFeatures_DISTANCE(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x) ;
	SVector& generateVertexFeatures_DISTANCE(unsigned i, const GraphClass& aG , int dag_h);
	void GenerateAbstractVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x);

};
#endif // CGRAPHPERCEPTRON_H
