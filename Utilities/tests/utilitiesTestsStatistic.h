// -*- C++ -*-
//
// utilitiesTestStatistic.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Utilities_Test_Statistic_H
#define HERWIG_Utilities_Test_Statistic_H

#include <boost/test/unit_test.hpp>

#include "Herwig/Utilities/Statistic.h"

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
    
BOOST_AUTO_TEST_SUITE_END()

#endif /* HERWIG_Utilities_Test_Statistic_H */