// -*- C++ -*-
//
// utilitiesTest.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "Herwig++/Utilities/Kinematics.h"
#include "Herwig++/Utilities/Maths.h"
#include "Herwig++/Utilities/Statistic.h"
#define BOOST_TEST_MODULE utilitiesTest
#include <boost/test/included/unit_test.hpp>


/*
 * Start of boost unit tests for Kinematics.h
 * 
 */
BOOST_AUTO_TEST_SUITE(utilitiesKinematicsTest)

/*
 * Boost unit tests
 * @todo solve problem with linking to ThePEG library
 */
/*
BOOST_AUTO_TEST_CASE(generateAngles)
{
  double flatMinusPiToPlusPi, flatNullToTwoPi;
  for(int i = 0; i < 100; ++i) {
    Herwig::Kinematics::generateAngles(flatMinusPiToPlusPi, flatNullToTwoPi);
    //BOOST_CHECK( -M_PI <= flatMinusPiToPlusPi );
    //BOOST_CHECK( flatMinusPiToPlusPi <= M_PI);
    //BOOST_CHECK( 0. <= flatNullToTwoPi);
    //BOOST_CHECK(flatNullToTwoPi <= M_PI);
  }
}
*/
    
BOOST_AUTO_TEST_SUITE_END()

//____________________________________________________________________________//

/*
 * Start of boost unit tests for Maths.h
 * 
 */
BOOST_AUTO_TEST_SUITE(utilitiesMathsTest)

/*
 * Boost unit tests
 * 
 */
BOOST_AUTO_TEST_CASE(dilogFunction)
{
  BOOST_CHECK_EQUAL(Herwig::Math::Li2(0.), 0.);
  BOOST_CHECK_CLOSE(Herwig::Math::Li2(-1).real(), -1./12. * M_PI * M_PI, 1e-5);
  BOOST_CHECK_CLOSE(Herwig::Math::Li2(-1).imag(), 0., 1e-5);
  BOOST_CHECK_CLOSE(Herwig::Math::Li2(1).real(), 1./6. * M_PI * M_PI, 1e-5);
  BOOST_CHECK_CLOSE(Herwig::Math::Li2(1).imag(), 0., 1e-5);
  
  BOOST_CHECK_CLOSE(Herwig::Math::Li2(Herwig::Complex(0.5, 1)).real(), 0.203354, 2e-4);
  BOOST_CHECK_CLOSE(Herwig::Math::Li2(Herwig::Complex(0.5, 1)).imag(), 1.13194, 1e-4);
  BOOST_CHECK_CLOSE(Herwig::Math::Li2(Herwig::Complex(-0.5, -1)).real(), -0.578503, 1e-4);
  BOOST_CHECK_CLOSE(Herwig::Math::Li2(Herwig::Complex(-0.5, -1)).imag(), -0.772714, 1e-4);
}

BOOST_AUTO_TEST_CASE(realDilogFunction)
{
  BOOST_CHECK_EQUAL(Herwig::Math::ReLi2(0.), 0.);
  BOOST_CHECK_CLOSE(Herwig::Math::ReLi2(-1), -1./12. * M_PI * M_PI, 1e-5);
  BOOST_CHECK_CLOSE(Herwig::Math::ReLi2(1), 1./6. * M_PI * M_PI, 1e-5);
}


BOOST_AUTO_TEST_CASE(angleZeroTo2Pi)
{
  BOOST_CHECK_EQUAL(Herwig::Math::angleZeroTo2Pi(0.), 0.);
  BOOST_CHECK_EQUAL(Herwig::Math::angleZeroTo2Pi(-0.5*M_PI), 1.5*M_PI);
  BOOST_CHECK_EQUAL(Herwig::Math::angleZeroTo2Pi(-2.5*M_PI), 1.5*M_PI);
  BOOST_CHECK_EQUAL(Herwig::Math::angleZeroTo2Pi( 2.5*M_PI), 0.5*M_PI);
  BOOST_CHECK(Herwig::Math::angleZeroTo2Pi( 2.5*M_PI) != 1*M_PI);
}

BOOST_AUTO_TEST_CASE(angleMinusPiToPi)
{
  BOOST_CHECK_EQUAL(Herwig::Math::angleMinusPiToPi(0.), 0.);
  BOOST_CHECK_EQUAL(Herwig::Math::angleMinusPiToPi(-0.5*M_PI), -0.5*M_PI);
  BOOST_CHECK_EQUAL(Herwig::Math::angleMinusPiToPi(-2.5*M_PI), -0.5*M_PI);
  BOOST_CHECK_EQUAL(Herwig::Math::angleMinusPiToPi( 2.5*M_PI), 0.5*M_PI);
  BOOST_CHECK(Herwig::Math::angleMinusPiToPi( 2.5*M_PI) != 1*M_PI);
}

BOOST_AUTO_TEST_CASE(median)
{
  std::vector<double> medianTest1;
  medianTest1.push_back(10);
  medianTest1.push_back(-1);
  medianTest1.push_back(5);
  BOOST_CHECK_EQUAL(Herwig::Math::median<double>(medianTest1), 5);
  
  std::vector<double> medianTest2;
  medianTest2.push_back(-10);
  medianTest2.push_back(-1);
  medianTest2.push_back(-5);
  medianTest2.push_back(-6);
  BOOST_CHECK_EQUAL(Herwig::Math::median<double>(medianTest2), -6);
}
    
BOOST_AUTO_TEST_SUITE_END()

//____________________________________________________________________________//

/*
 * Fixture which defines the common variables for testing Statistic class
 */
struct FixStatistic1 {
  FixStatistic1() : statisticDefault(), statisticTest()  
  {BOOST_TEST_MESSAGE( "setup fixture for utilitiesStatisticTest" ); }
  
  ~FixStatistic1()  { BOOST_TEST_MESSAGE( "teardown fixture for utilitiesStatisticTest" ); }

  Herwig::Statistic statisticDefault;
  Herwig::Statistic statisticTest;
};

/*
 * Start of boost unit tests for Statistic.h
 * 
 */
BOOST_FIXTURE_TEST_SUITE(utilitiesStatisticTest, FixStatistic1 )

/*
 * Boost unit tests 
 * 
 */
BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  BOOST_CHECK_EQUAL(statisticDefault.numberOfPoints(), static_cast<unsigned int>(0));
  BOOST_CHECK_EQUAL(statisticDefault.total(), 0.);
  BOOST_CHECK_EQUAL(statisticDefault.mean(), 0.);
  BOOST_CHECK_EQUAL(statisticDefault.minimum(), -1e100);
  BOOST_CHECK_EQUAL(statisticDefault.maximum(), 1e100);
  BOOST_CHECK_EQUAL(statisticDefault.var(), 0);
  BOOST_CHECK_EQUAL(statisticDefault.mean_var(), 0);
}

BOOST_AUTO_TEST_CASE(operations)
{
  statisticTest += 2;
  BOOST_CHECK_EQUAL(statisticTest.minimum(), 2);
  BOOST_CHECK_EQUAL(statisticTest.maximum(), 2);
  
  statisticTest += -2;
  BOOST_CHECK_EQUAL(statisticTest.numberOfPoints(), static_cast<unsigned int>(2));
  BOOST_CHECK_EQUAL(statisticTest.total(), 0.);
  BOOST_CHECK_EQUAL(statisticTest.mean(), 0.);
  BOOST_CHECK_EQUAL(statisticTest.minimum(), -2);
  BOOST_CHECK_EQUAL(statisticTest.maximum(), 2);
  BOOST_CHECK_EQUAL(statisticTest.var(), 8);
  BOOST_CHECK_EQUAL(statisticTest.stdDev(), sqrt(8));
  BOOST_CHECK_EQUAL(statisticTest.mean_var(), 4);
  BOOST_CHECK_EQUAL(statisticTest.mean_stdDev(), 2);
  
  statisticTest += 2;
  BOOST_CHECK_EQUAL(statisticTest.numberOfPoints(), static_cast<unsigned int>(3));
  BOOST_CHECK_EQUAL(statisticTest.total(), 2.);
  BOOST_CHECK_EQUAL(statisticTest.mean(), 2./3.);
  BOOST_CHECK_EQUAL(statisticTest.minimum(), -2);
  BOOST_CHECK_EQUAL(statisticTest.maximum(), 2);
  BOOST_CHECK_EQUAL(statisticTest.var(), 16./3.);
  BOOST_CHECK_EQUAL(statisticTest.stdDev(), sqrt(16./3.));
  BOOST_CHECK_EQUAL(statisticTest.mean_var(), 16./9.);
  BOOST_CHECK_EQUAL(statisticTest.mean_stdDev(), sqrt(16./9.));
  
  statisticTest += 4.5;
  BOOST_CHECK_EQUAL(statisticTest.minimum(), -2);
  BOOST_CHECK_EQUAL(statisticTest.maximum(), 4.5);
  
  statisticTest += -3.5;
  BOOST_CHECK_EQUAL(statisticTest.minimum(), -3.5);
  BOOST_CHECK_EQUAL(statisticTest.maximum(), 4.5);   
}    

BOOST_AUTO_TEST_CASE(fail)
{
  BOOST_FAIL("Ende");
}

    
BOOST_AUTO_TEST_SUITE_END()

//____________________________________________________________________________//

