// -*- C++ -*-
//
// Herwig++.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "HerwigRun.h"
#include <iostream>
#include "ThePEG/Utilities/Exception.h"

int main(int argc, char * argv[]) {
  // HerwigRun's constructor does all the work
  try {
    Herwig::HerwigRun hw(argc, argv);
    return hw.good() ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  catch ( ThePEG::Exception & e ) {
    return EXIT_FAILURE;
  }
  catch ( std::exception & e ) {
    std::cerr << argv[0] << ": " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  catch (...) {
    std::cerr << argv[0] << ": Unknown exception caught.\n";
    return EXIT_FAILURE;
  }
  
}
