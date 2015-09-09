// -*- C++ -*-
//
// utilitiesTestGlobalFixture.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include <boost/test/unit_test.hpp>

#include "ThePEG/Repository/StandardRandom.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Config/Unitsystem.h"

struct FixGlobal1 {
  FixGlobal1() {
    BOOST_TEST_MESSAGE( "setup global fixture for utilitiesTest" ); 
    
    // Initialize randomNumberGenerator
    ThePEG::StandardRandom* randomNumberStandardGenerator = new ThePEG::StandardRandom();
    new ThePEG::UseRandom(randomNumberStandardGenerator);
  }
  
  ~FixGlobal1()  { BOOST_TEST_MESSAGE( "teardown global fixture for utilitiesTest" ); }
};

BOOST_GLOBAL_FIXTURE(FixGlobal1);
