// -*- C++ -*-
//
// utilitiesTestMaths.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Utilities_Test_Maths_H
#define HERWIG_Utilities_Test_Maths_H

#include <boost/test/unit_test.hpp>

#include "Herwig/Utilities/Maths.h"

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

#endif /* HERWIG_Utilities_Test_Maths_H */