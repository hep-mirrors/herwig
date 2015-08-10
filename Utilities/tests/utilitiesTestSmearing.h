// -*- C++ -*-
//
// utilitiesTestSmearing.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Utilities_Test_Smearing_H
#define HERWIG_Utilities_Test_Smearing_H

#include <boost/test/unit_test.hpp>

#include "Herwig++/Utilities/Smearing.h"
#include <iostream>
#include <fstream>

/*
 * Start of boost unit tests for Smearing.h
 * 
 */
BOOST_AUTO_TEST_SUITE(utilitiesSmearingTest)

/*
 * Boost unit tests 
 * 
 */
BOOST_AUTO_TEST_CASE(azimuthalSmearing)
{
  double vectorXvalue, vectorYvalue;
  double vectorLengthCalculated;
  unsigned int numberOfTrials = 10000;
  for(unsigned int vectorLength = 1; vectorLength < 4; ++vectorLength) {
    // Variables for statistical analysis
    double meanVectorLength = 0, meanVectorX = 0, meanVectorY = 0;
    unsigned int inverseSamplingEfficiency = 0;
    for(unsigned int i = 0; i < numberOfTrials; ++i) {
      if(Herwig::Smearing::azimuthalSmearing(vectorLength, vectorXvalue, vectorYvalue)) {
	vectorLengthCalculated = sqrt(vectorXvalue*vectorXvalue + vectorYvalue*vectorYvalue);
	BOOST_CHECK_CLOSE(vectorLengthCalculated, vectorLength, 1e-06);
	meanVectorLength += vectorLengthCalculated;
	meanVectorX += vectorXvalue;
	meanVectorY += vectorYvalue;
      }
      else {
	inverseSamplingEfficiency++;
      }
    }
  // Normalize mean values
  meanVectorLength = meanVectorLength/static_cast<double>(numberOfTrials - inverseSamplingEfficiency);
  meanVectorX = meanVectorX/static_cast<double>(numberOfTrials - inverseSamplingEfficiency);
  meanVectorY = meanVectorY/static_cast<double>(numberOfTrials - inverseSamplingEfficiency);
  
  BOOST_CHECK_CLOSE(meanVectorLength , vectorLength, 0.001);
  BOOST_CHECK(meanVectorX < 1e-1 && meanVectorX > -1e-1);
  BOOST_CHECK(meanVectorY < 1e-1 && meanVectorY > -1e-1);
  BOOST_CHECK(inverseSamplingEfficiency < 0.25 * numberOfTrials);
  }
}

BOOST_AUTO_TEST_CASE(gaussianSmearing)
{
  std::ofstream myfile;
  myfile.open ("out2.txt");
  
  double gaussianValue;
  unsigned int numberOfTrials = 1000000;
  for(unsigned int gaussianMean = 1; gaussianMean < 2; ++gaussianMean) {
    for(double gaussianSigma = 0.25; gaussianSigma < 0.35; gaussianSigma += 0.1) {
      // Variables for statistical analysis
      double meanValue = 0, meanSigma = 0;
      unsigned int inverseSamplingEfficiency = 0;
      for(unsigned int i = 0; i < numberOfTrials; ++i) {
	if(Herwig::Smearing::gaussianSmearing(gaussianMean, gaussianSigma, gaussianValue)) {
	  // Check for 8 sigma deviation
	  BOOST_CHECK(gaussianValue < (gaussianMean + 8 * gaussianSigma));
	  //std::cout << gaussianValue << " " << gaussianMean << " " << gaussianSigma << " " << (gaussianMean + 8 * gaussianSigma) << std::endl;
	  myfile << gaussianValue << std::endl;
	  BOOST_CHECK(gaussianValue > (gaussianMean - 8 * gaussianSigma));
	  meanValue += gaussianValue;
	  meanSigma += std::abs(gaussianValue - gaussianMean);
	}
	else {
	  inverseSamplingEfficiency++;
	}
      }
    // Normalize mean values
    meanValue = meanValue/static_cast<double>(numberOfTrials - inverseSamplingEfficiency);
    meanSigma = meanSigma/static_cast<double>(numberOfTrials - inverseSamplingEfficiency);
    std::cout << meanValue << " " << meanSigma << std::endl;
    BOOST_CHECK_CLOSE(meanValue , gaussianMean, 0.5);
    BOOST_CHECK_CLOSE(meanSigma , gaussianSigma, 0.5);
    BOOST_CHECK(inverseSamplingEfficiency < 0.25 * numberOfTrials);
    }
  }
}



BOOST_AUTO_TEST_SUITE_END()

#endif /* HERWIG_Utilities_Test_Smearing_H */