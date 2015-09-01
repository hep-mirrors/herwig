// -*- C++ -*-
//
// utilitiesTestKinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Utilities_Test_Kinematics_H
#define HERWIG_Utilities_Test_Kinematics_H

#include <boost/test/unit_test.hpp>

#include "Herwig/Utilities/Kinematics.h"

#include "ThePEG/Config/Unitsystem.h"


struct FixKinematics1 {
  FixKinematics1()  
  {BOOST_TEST_MESSAGE( "setup fixture for utilitiesKinematicsTestTest" );  }
  
  ~FixKinematics1()  { BOOST_TEST_MESSAGE( "teardown fixture for utilitiesKinematicsTest" ); }
};

/*
 * Start of boost unit tests for Kinematics.h
 *
 * @todo Implement unit test for threeBodyDecay
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
    
BOOST_AUTO_TEST_CASE(unitDirection)
{
  BOOST_CHECK_EQUAL( Herwig::Kinematics::unitDirection(1.1, -1), ThePEG::Units::Axis());
  BOOST_CHECK_EQUAL( Herwig::Kinematics::unitDirection(-1.1, -1), ThePEG::Units::Axis());
  BOOST_CHECK_EQUAL( Herwig::Kinematics::unitDirection(1.1, 0), ThePEG::Units::Axis());
  BOOST_CHECK_EQUAL( Herwig::Kinematics::unitDirection(-1.1, -1), ThePEG::Units::Axis());
  
  BOOST_CHECK_EQUAL( Herwig::Kinematics::unitDirection(1, 0), ThePEG::Units::Axis(0, 0, 1));
  BOOST_CHECK_EQUAL( Herwig::Kinematics::unitDirection(1, M_PI/2.), ThePEG::Units::Axis(0, 0, 1));
  
  BOOST_CHECK(Herwig::Kinematics::unitDirection(0, M_PI/2).almostEqual(ThePEG::Units::Axis(0, 1, 0), 0.001));
  BOOST_CHECK_EQUAL( Herwig::Kinematics::unitDirection(0, 0), ThePEG::Units::Axis(1, 0, 0));
}
    
BOOST_AUTO_TEST_CASE(pstarTwoBodyDecay)
{
  BOOST_CHECK_EQUAL(Herwig::Kinematics::pstarTwoBodyDecay(ThePEG::Units::Energy(-100*ThePEG::Units::GeV), ThePEG::Units::Energy(60*ThePEG::Units::GeV), ThePEG::Units::Energy(60*ThePEG::Units::GeV))/ThePEG::Units::GeV, ThePEG::Units::Energy(0*ThePEG::Units::GeV)/ThePEG::Units::GeV);
  BOOST_CHECK_EQUAL(Herwig::Kinematics::pstarTwoBodyDecay(ThePEG::Units::Energy(100*ThePEG::Units::GeV), ThePEG::Units::Energy(-40*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV))/ThePEG::Units::GeV, ThePEG::Units::Energy(0*ThePEG::Units::GeV)/ThePEG::Units::GeV);
  BOOST_CHECK_EQUAL(Herwig::Kinematics::pstarTwoBodyDecay(ThePEG::Units::Energy(100*ThePEG::Units::GeV), ThePEG::Units::Energy(-40*ThePEG::Units::GeV), ThePEG::Units::Energy(-40*ThePEG::Units::GeV))/ThePEG::Units::GeV, ThePEG::Units::Energy(0*ThePEG::Units::GeV)/ThePEG::Units::GeV);
  
  BOOST_CHECK_EQUAL(Herwig::Kinematics::pstarTwoBodyDecay(ThePEG::Units::Energy(100*ThePEG::Units::GeV), ThePEG::Units::Energy(60*ThePEG::Units::GeV), ThePEG::Units::Energy(60*ThePEG::Units::GeV))/ThePEG::Units::GeV, ThePEG::Units::Energy(0*ThePEG::Units::GeV)/ThePEG::Units::GeV);
  BOOST_CHECK_EQUAL(Herwig::Kinematics::pstarTwoBodyDecay(ThePEG::Units::Energy(100*ThePEG::Units::GeV), ThePEG::Units::Energy(50*ThePEG::Units::GeV), ThePEG::Units::Energy(50*ThePEG::Units::GeV))/ThePEG::Units::GeV, ThePEG::Units::Energy(0*ThePEG::Units::GeV)/ThePEG::Units::GeV);
  
  BOOST_CHECK_EQUAL(Herwig::Kinematics::pstarTwoBodyDecay(ThePEG::Units::Energy(10*ThePEG::Units::GeV), ThePEG::Units::Energy(6*ThePEG::Units::GeV), ThePEG::Units::Energy(3*ThePEG::Units::GeV))/ThePEG::Units::GeV, ThePEG::Units::Energy(std::sqrt(19*91)/20.*ThePEG::Units::GeV)/ThePEG::Units::GeV);
  BOOST_CHECK_EQUAL(Herwig::Kinematics::pstarTwoBodyDecay(ThePEG::Units::Energy(10*ThePEG::Units::GeV), ThePEG::Units::Energy(3*ThePEG::Units::GeV), ThePEG::Units::Energy(6*ThePEG::Units::GeV))/ThePEG::Units::GeV, ThePEG::Units::Energy(std::sqrt(19*91)/20.*ThePEG::Units::GeV)/ThePEG::Units::GeV);
}
    
BOOST_AUTO_TEST_CASE(twoBodyDecay1)
{
 ThePEG::Units::Lorentz5Momentum decayProductOne(ThePEG::Units::GeV);
 ThePEG::Units::Lorentz5Momentum decayProductTwo(ThePEG::Units::GeV);
 BOOST_CHECK(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), ThePEG::Units::Axis(), decayProductOne, decayProductTwo));
 BOOST_CHECK(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(50*ThePEG::Units::GeV), ThePEG::Units::Energy(50*ThePEG::Units::GeV), ThePEG::Units::Axis(), decayProductOne, decayProductTwo));
 
 BOOST_CHECK(!(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(60*ThePEG::Units::GeV), ThePEG::Units::Energy(60*ThePEG::Units::GeV), ThePEG::Units::Axis(), decayProductOne, decayProductTwo))); 
 BOOST_CHECK(!(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(-100*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), ThePEG::Units::Axis(), decayProductOne, decayProductTwo)));
 BOOST_CHECK(!(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(-40*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), ThePEG::Units::Axis(), decayProductOne, decayProductTwo)));
 BOOST_CHECK(!(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), ThePEG::Units::Energy(-40*ThePEG::Units::GeV), ThePEG::Units::Axis(), decayProductOne, decayProductTwo)));
 
 Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(50*ThePEG::Units::GeV), ThePEG::Units::Energy(50*ThePEG::Units::GeV), ThePEG::Units::Axis(1,0,0), decayProductOne, decayProductTwo);
 BOOST_CHECK_EQUAL(decayProductOne/ThePEG::Units::GeV, ThePEG::Units::Lorentz5Momentum(50*ThePEG::Units::GeV)/ThePEG::Units::GeV);
 BOOST_CHECK_EQUAL(decayProductTwo/ThePEG::Units::GeV, ThePEG::Units::Lorentz5Momentum(50*ThePEG::Units::GeV)/ThePEG::Units::GeV); 
 
 Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(10*ThePEG::Units::GeV), ThePEG::Units::Energy(6*ThePEG::Units::GeV), ThePEG::Units::Energy(3*ThePEG::Units::GeV), ThePEG::Units::Axis(1,0,0), decayProductOne, decayProductTwo);
 BOOST_CHECK_EQUAL(decayProductOne/ThePEG::Units::GeV, ThePEG::Units::Lorentz5Momentum(6*ThePEG::Units::GeV, ThePEG::Units::Momentum3(std::sqrt(19*91)/20.*ThePEG::Units::GeV, ThePEG::ZERO, ThePEG::ZERO))/ThePEG::Units::GeV);
 BOOST_CHECK_EQUAL(decayProductTwo/ThePEG::Units::GeV, ThePEG::Units::Lorentz5Momentum(3*ThePEG::Units::GeV, ThePEG::Units::Momentum3(-(std::sqrt(19*91)/20.*ThePEG::Units::GeV), ThePEG::ZERO, ThePEG::ZERO))/ThePEG::Units::GeV);
}

BOOST_AUTO_TEST_CASE(twoBodyDecay2)
{
 ThePEG::Units::Lorentz5Momentum decayProductOne(ThePEG::Units::GeV);
 ThePEG::Units::Lorentz5Momentum decayProductTwo(ThePEG::Units::GeV);
 BOOST_CHECK(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), 1, M_PI/2., decayProductOne, decayProductTwo));
 BOOST_CHECK(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(50*ThePEG::Units::GeV), ThePEG::Units::Energy(50*ThePEG::Units::GeV), 1, M_PI/2., decayProductOne, decayProductTwo));
 
 BOOST_CHECK(!(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(60*ThePEG::Units::GeV), ThePEG::Units::Energy(60*ThePEG::Units::GeV), 1, M_PI/2., decayProductOne, decayProductTwo))); 
 BOOST_CHECK(!(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(-100*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), 1, M_PI/2., decayProductOne, decayProductTwo)));
 BOOST_CHECK(!(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(-40*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), 1, M_PI/2., decayProductOne, decayProductTwo)));
 BOOST_CHECK(!(Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(40*ThePEG::Units::GeV), ThePEG::Units::Energy(-40*ThePEG::Units::GeV), 1, M_PI/2., decayProductOne, decayProductTwo)));
 
 Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(100*ThePEG::Units::GeV), ThePEG::Units::Energy(50*ThePEG::Units::GeV), ThePEG::Units::Energy(50*ThePEG::Units::GeV), 1, M_PI/2., decayProductOne, decayProductTwo);
 BOOST_CHECK_EQUAL(decayProductOne/ThePEG::Units::GeV, ThePEG::Units::Lorentz5Momentum(50*ThePEG::Units::GeV)/ThePEG::Units::GeV);
 BOOST_CHECK_EQUAL(decayProductTwo/ThePEG::Units::GeV, ThePEG::Units::Lorentz5Momentum(50*ThePEG::Units::GeV)/ThePEG::Units::GeV); 
 
 Herwig::Kinematics::twoBodyDecay(ThePEG::Units::Lorentz5Momentum(10*ThePEG::Units::GeV), ThePEG::Units::Energy(6*ThePEG::Units::GeV), ThePEG::Units::Energy(3*ThePEG::Units::GeV), 1, M_PI/2., decayProductOne, decayProductTwo);
 BOOST_CHECK_EQUAL(decayProductOne/ThePEG::Units::GeV, ThePEG::Units::Lorentz5Momentum(6*ThePEG::Units::GeV, ThePEG::Units::Momentum3(ThePEG::ZERO, ThePEG::ZERO,   std::sqrt(19*91)/20.*ThePEG::Units::GeV))/ThePEG::Units::GeV);
 BOOST_CHECK_EQUAL(decayProductTwo/ThePEG::Units::GeV, ThePEG::Units::Lorentz5Momentum(3*ThePEG::Units::GeV, ThePEG::Units::Momentum3(ThePEG::ZERO, ThePEG::ZERO, -(std::sqrt(19*91)/20.*ThePEG::Units::GeV)))/ThePEG::Units::GeV);
}


BOOST_AUTO_TEST_SUITE_END()

#endif /* HERWIG_Utilities_Test_Kinematics_H */
