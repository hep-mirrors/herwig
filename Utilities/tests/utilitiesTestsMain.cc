// -*- C++ -*-
//
// utilitiesTestMain.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

/**
 * The following part should be included only once. 
*/
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE utilitiesTest

/**
 * Include global fixture
 * 
 * Global fixture initializes the randomNumber generator
 */
#include "Herwig/Utilities/tests/utilitiesTestsGlobalFixture.h"

/**
 * Include here the sub tests
 */
#include "Herwig/Utilities/tests/utilitiesTestsKinematics.h"
#include "Herwig/Utilities/tests/utilitiesTestMaths.h"
#include "Herwig/Utilities/tests/utilitiesTestsStatistic.h"


/**
 * Debug and development part
 */
BOOST_AUTO_TEST_CASE(fail)
{
  //BOOST_FAIL("Ende");
}