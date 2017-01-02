/* -*- mode:c++ -*-
 * FeatureGenerator.h
 *
 *  Created on: Jul 23, 2010
 *      Author: Paolo Frasconi
 */

#ifndef FEATUREGENERATOR_H_
#define FEATUREGENERATOR_H_

#include <string>
#include <iostream>
#include <sstream>
#include "vectors.h"

class FeatureGenerator {
protected:
  std::string feature_generator_id;

public:
  FeatureGenerator(std::string feature_generator_id = "void_feature_generator") {
    this->feature_generator_id = feature_generator_id;
  };
  virtual ~FeatureGenerator() {};
  virtual std::string str(void) {
    std::stringstream oss;
    oss << "FeatureGenerator: >" << feature_generator_id << "<" << std::endl;
    return oss.str();
  }
  virtual void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList=vector<unsigned>()) = 0;
};

#endif /* FEATUREGENERATOR_H_ */

