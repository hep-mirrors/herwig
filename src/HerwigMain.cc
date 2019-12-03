// -*- C++ -*-
//
// HerwigMain.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "HerwigCLI.h"
#include "Herwig/API/HerwigAPI.h"
#include <iostream>
#include <cstdlib>
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Utilities/ColourOutput.h"



int main(int argc, char * argv[]) {
 
  try {

    // read in command line options
    const Herwig::HerwigCLI cl = Herwig::HerwigCLI(argc, argv);
  
    // Call program switches according to runMode
    switch ( cl.runMode() ) {
    case Herwig::RunMode::INIT:        Herwig::API::init(cl);       break;
    case Herwig::RunMode::READ:        Herwig::API::read(cl);       break;
    case Herwig::RunMode::BUILD:       Herwig::API::build(cl);      break;
    case Herwig::RunMode::INTEGRATE:   Herwig::API::integrate(cl);  break;
    case Herwig::RunMode::MERGEGRIDS:  Herwig::API::mergegrids(cl); break;
    case Herwig::RunMode::RUN:         Herwig::API::run(cl);        break;
    case Herwig::RunMode::ERROR:       
      std::cerr << ThePEG::ANSI::red 
    	        << "Error during read in of command line parameters.\n"
                << "Program execution will stop now."
    	        << ThePEG::ANSI::reset;
      return EXIT_FAILURE;
    default:          		     
      cl.quitWithHelp();
    }

    return EXIT_SUCCESS;

  }
  catch ( ThePEG::Exception & e ) {
    std::cerr << ThePEG::ANSI::red 
    	      << argv[0] << ": ThePEG::Exception caught.\n"
              << e.what() << '\n'
              << ThePEG::ANSI::yellow
      	      << "See logfile for details.\n"
      	      << ThePEG::ANSI::reset;
    return EXIT_FAILURE;
  }
  catch ( std::exception & e ) {
    std::cerr << ThePEG::ANSI::red 
    	      << argv[0] << ": " << e.what() << '\n'
    	      << ThePEG::ANSI::reset;
    return EXIT_FAILURE;
  }
  catch ( const char* what ) {
    std::cerr << ThePEG::ANSI::red 
    	      << argv[0] << ": caught exception: "
	      << what << "\n"
    	      << ThePEG::ANSI::reset;
    return EXIT_FAILURE;
  }
  catch ( ... ) {
    std::cerr << ThePEG::ANSI::red 
    	      << argv[0] << ": Unknown exception caught.\n"
    	      << ThePEG::ANSI::reset;
    return EXIT_FAILURE;
  }  
} 
