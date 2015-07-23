// -*- C++ -*-
//
// utilitiesTestKinematics.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Utilities_Test_Kinematics_H
#define HERWIG_Utilities_Test_Kinematics_H

#include <boost/test/unit_test.hpp>

#include "Herwig++/Utilities/Kinematics.h"

struct FixKinematics1 {
  FixKinematics1()  
  {BOOST_TEST_MESSAGE( "setup fixture for utilitiesKinematicsTestTest" );  }
  
  ~FixKinematics1()  { BOOST_TEST_MESSAGE( "teardown fixture for utilitiesKinematicsTest" ); }
};

/*
 * Start of boost unit tests for Kinematics.h
 * 
 */
BOOST_AUTO_TEST_SUITE(utilitiesKinematicsTest)

/*
 * Boost unit tests
 *
 */
BOOST_AUTO_TEST_CASE(generateAnglestest)
{
  double flatMinusPiToPlusPi, flatNullToTwoPi;
  for(int i = 0; i < 100; ++i) {
    Herwig::Kinematics::generateAngles(flatMinusPiToPlusPi, flatNullToTwoPi);
    BOOST_CHECK( -M_PI <= flatMinusPiToPlusPi );
    BOOST_CHECK( flatMinusPiToPlusPi <= M_PI);
    BOOST_CHECK( 0. <= flatNullToTwoPi);
    BOOST_CHECK(flatNullToTwoPi <= 2*M_PI);
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* HERWIG_Utilities_Test_Kinematics_H */