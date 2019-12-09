// -*- C++ -*-
//
// utilitiesTestKinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Utilities_Test_Kinematics_H
#define HERWIG_Utilities_Test_Kinematics_H

#include <boost/test/unit_test.hpp>

#include "Herwig/Utilities/Kinematics.h"

#include "ThePEG/Config/Unitsystem.h"

using namespace Herwig::Kinematics;
using namespace ThePEG::Units;

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
BOOST_AUTO_TEST_CASE(generateAnglesTest)
{
  double flatMinusPiToPlusPi, flatNullToTwoPi;
  for(int i = 0; i < 100; ++i) {
    generateAngles(flatMinusPiToPlusPi, flatNullToTwoPi);
    BOOST_CHECK( -M_PI <= flatMinusPiToPlusPi );
    BOOST_CHECK( flatMinusPiToPlusPi <= M_PI);
    BOOST_CHECK( 0. <= flatNullToTwoPi);
    BOOST_CHECK(flatNullToTwoPi <= 2*M_PI);
  }
}
    
BOOST_AUTO_TEST_CASE(unitDirectionTest)
{
  using namespace Herwig::Kinematics;
  BOOST_CHECK_EQUAL( unitDirection(1.1, -1), Axis() );
  BOOST_CHECK_EQUAL( unitDirection(-1.1, -1), Axis() );
  BOOST_CHECK_EQUAL( unitDirection(1.1, 0), Axis() );
  BOOST_CHECK_EQUAL( unitDirection(-1.1, -1), Axis() );
  
  BOOST_CHECK_EQUAL( unitDirection(1, 0), Axis(0, 0, 1) );
  BOOST_CHECK_EQUAL( unitDirection(1, M_PI/2.), Axis(0, 0, 1) );
  
  BOOST_CHECK(unitDirection(0, M_PI/2).almostEqual(Axis(0, 1, 0), 0.001) );
  BOOST_CHECK_EQUAL( unitDirection(0, 0), Axis(1, 0, 0) );
}
    
BOOST_AUTO_TEST_CASE(pstarTwoBodyDecayTest)
{
  BOOST_CHECK_EQUAL(pstarTwoBodyDecay(Energy(-100*GeV), Energy(60*GeV), Energy(60*GeV))/GeV, Energy(0*GeV)/GeV);
  BOOST_CHECK_EQUAL(pstarTwoBodyDecay(Energy(100*GeV), Energy(-40*GeV), Energy(40*GeV))/GeV, Energy(0*GeV)/GeV);
  BOOST_CHECK_EQUAL(pstarTwoBodyDecay(Energy(100*GeV), Energy(-40*GeV), Energy(-40*GeV))/GeV, Energy(0*GeV)/GeV);
  
  BOOST_CHECK_EQUAL(pstarTwoBodyDecay(Energy(100*GeV), Energy(60*GeV), Energy(60*GeV))/GeV, Energy(0*GeV)/GeV);
  BOOST_CHECK_EQUAL(pstarTwoBodyDecay(Energy(100*GeV), Energy(50*GeV), Energy(50*GeV))/GeV, Energy(0*GeV)/GeV);
  
  BOOST_CHECK_EQUAL(pstarTwoBodyDecay(Energy(10*GeV), Energy(6*GeV), Energy(3*GeV))/GeV, Energy(std::sqrt(19*91)/20.*GeV)/GeV);
  BOOST_CHECK_EQUAL(pstarTwoBodyDecay(Energy(10*GeV), Energy(3*GeV), Energy(6*GeV))/GeV, Energy(std::sqrt(19*91)/20.*GeV)/GeV);
}
    
BOOST_AUTO_TEST_CASE(twoBodyDecayTest1)
{
 Lorentz5Momentum decayProductOne(GeV);
 Lorentz5Momentum decayProductTwo(GeV);
 BOOST_CHECK(twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(40*GeV), Energy(40*GeV), Axis(), decayProductOne, decayProductTwo));
 BOOST_CHECK(twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(50*GeV), Energy(50*GeV), Axis(), decayProductOne, decayProductTwo));
 
 BOOST_CHECK(!(twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(60*GeV), Energy(60*GeV), Axis(), decayProductOne, decayProductTwo))); 
 BOOST_CHECK(!(twoBodyDecay(Lorentz5Momentum(-100*GeV), Energy(40*GeV), Energy(40*GeV), Axis(), decayProductOne, decayProductTwo)));
 BOOST_CHECK(!(twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(-40*GeV), Energy(40*GeV), Axis(), decayProductOne, decayProductTwo)));
 BOOST_CHECK(!(twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(40*GeV), Energy(-40*GeV), Axis(), decayProductOne, decayProductTwo)));
 
 twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(50*GeV), Energy(50*GeV), Axis(1,0,0), decayProductOne, decayProductTwo);
 BOOST_CHECK_EQUAL(decayProductOne/GeV, Lorentz5Momentum(50*GeV)/GeV);
 BOOST_CHECK_EQUAL(decayProductTwo/GeV, Lorentz5Momentum(50*GeV)/GeV); 
 
 twoBodyDecay(Lorentz5Momentum(10*GeV), Energy(6*GeV), Energy(3*GeV), Axis(1,0,0), decayProductOne, decayProductTwo);
 BOOST_CHECK_EQUAL(decayProductOne/GeV, Lorentz5Momentum(6*GeV, Momentum3(std::sqrt(19*91)/20.*GeV, ThePEG::ZERO, ThePEG::ZERO))/GeV);
 BOOST_CHECK_EQUAL(decayProductTwo/GeV, Lorentz5Momentum(3*GeV, Momentum3(-(std::sqrt(19*91)/20.*GeV), ThePEG::ZERO, ThePEG::ZERO))/GeV);
}

BOOST_AUTO_TEST_CASE(twoBodyDecayTest2)
{
 Lorentz5Momentum decayProductOne(GeV);
 Lorentz5Momentum decayProductTwo(GeV);
 BOOST_CHECK(twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(40*GeV), Energy(40*GeV), 1, M_PI/2., decayProductOne, decayProductTwo));
 BOOST_CHECK(twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(50*GeV), Energy(50*GeV), 1, M_PI/2., decayProductOne, decayProductTwo));
 
 BOOST_CHECK(!(twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(60*GeV), Energy(60*GeV), 1, M_PI/2., decayProductOne, decayProductTwo))); 
 BOOST_CHECK(!(twoBodyDecay(Lorentz5Momentum(-100*GeV), Energy(40*GeV), Energy(40*GeV), 1, M_PI/2., decayProductOne, decayProductTwo)));
 BOOST_CHECK(!(twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(-40*GeV), Energy(40*GeV), 1, M_PI/2., decayProductOne, decayProductTwo)));
 BOOST_CHECK(!(twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(40*GeV), Energy(-40*GeV), 1, M_PI/2., decayProductOne, decayProductTwo)));
 
 twoBodyDecay(Lorentz5Momentum(100*GeV), Energy(50*GeV), Energy(50*GeV), 1, M_PI/2., decayProductOne, decayProductTwo);
 BOOST_CHECK_EQUAL(decayProductOne/GeV, Lorentz5Momentum(50*GeV)/GeV);
 BOOST_CHECK_EQUAL(decayProductTwo/GeV, Lorentz5Momentum(50*GeV)/GeV); 
 
 twoBodyDecay(Lorentz5Momentum(10*GeV), Energy(6*GeV), Energy(3*GeV), 1, M_PI/2., decayProductOne, decayProductTwo);
 BOOST_CHECK_EQUAL(decayProductOne/GeV, Lorentz5Momentum(6*GeV, Momentum3(ThePEG::ZERO, ThePEG::ZERO,   std::sqrt(19*91)/20.*GeV))/GeV);
 BOOST_CHECK_EQUAL(decayProductTwo/GeV, Lorentz5Momentum(3*GeV, Momentum3(ThePEG::ZERO, ThePEG::ZERO, -(std::sqrt(19*91)/20.*GeV)))/GeV);
}


BOOST_AUTO_TEST_SUITE_END()

#endif /* HERWIG_Utilities_Test_Kinematics_H */
