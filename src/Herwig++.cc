// -*- C++ -*-
//
// Herwig++.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "herwigopts.h"
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Utilities/DynamicLoader.h>
#include <ThePEG/Repository/Repository.h>
#include <ThePEG/Utilities/Exception.h>
#include <ThePEG/Utilities/Debug.h>
#include <iostream>

using namespace ThePEG;

void printUsageAndExit();

void HerwigInit(string infile, string reponame);
void HerwigRead(string reponame, string runname);
void HerwigRun(string runname, int seed, long N, 
	       bool tics, bool resume, bool keepid);



int main(int argc, char * argv[]) {
 
  try {

    // read command line options 
    gengetopt_args_info args_info;
    if ( cmdline_parser( argc, argv, &args_info ) != 0 ) {
      std::cerr << "Could not parse command line.\n";
      return EXIT_FAILURE;
    }
 
    // require one command
    if ( args_info.inputs_num < 1 )
      printUsageAndExit();

    // Interpret command status
    enum { INIT, READ, RUN } status;
    std::string runType = args_info.inputs[0];
    if ( runType == "init" )
      status = INIT;
    else if ( runType == "read" ) 
      status = READ;
    else if ( runType == "run" )  
      status = RUN;
    else                          
      printUsageAndExit();

    // Use second argument as input- or runfile name
    string runname;
    if ( args_info.inputs_num > 1 )
      runname = args_info.inputs[1];

    // If status is RUN, we need a runname
    if ( status == RUN && runname.empty() ) {
      cerr << "Error: You need to supply a runfile name.\n";
      printUsageAndExit();
    }

    // Defaults for these filenames are set in the ggo file
    std::string reponame = args_info.repo_arg;
    std::string infile = args_info.init_arg;

    // Number of events
    long N = -1;
    if ( args_info.numevents_given )
      N = args_info.numevents_arg;

    // RNG seed
    int seed = 0;
    if ( args_info.seed_given )
      seed = args_info.seed_arg;

    // Library search path for dlopen()
    for ( size_t i = 0; i < args_info.append_given; ++i )
      DynamicLoader::appendPath( args_info.append_arg[i] );
    for ( size_t i = 0; i < args_info.prepend_given; ++i )
      DynamicLoader::prependPath( args_info.prepend_arg[i] );

    // Debugging level
    if ( args_info.debug_given )
      Debug::setDebug( args_info.debug_arg );
    if ( Debug::level ) 
      Debug::unmaskFpuErrors();

    // Exit-on-error flag
    if ( args_info.exitonerror_flag )
      Repository::exitOnError() = 1;

    // Tics
    bool tics = true;
    if ( args_info.quiet_flag )
      tics = false;

    // Resume
    bool resume = false;
    if ( args_info.resume_flag )
      resume = true;

    // Keep id
    bool keepid = false;
    if ( args_info.keepid_flag )
      keepid = true;

    // *** End of command line parsing ***
   
    // Call mode
    switch ( status ) {
    case INIT:  HerwigInit( infile, reponame ); break;
    case READ:  HerwigRead( reponame, runname ); break;
    case RUN:   HerwigRun( runname, seed, N, 
			   tics, resume, keepid );  break;
    default:    printUsageAndExit();
    }


    Repository::cleanup(); 
    cmdline_parser_free( &args_info );
    return EXIT_SUCCESS;

  }
  catch ( ThePEG::Exception & e ) {
    std::cerr << argv[0] << ": ThePEG::Exception caught. "
	      << "See logfile for details.\n";
    Repository::cleanup(); 
    return EXIT_FAILURE;
  }
  catch ( std::exception & e ) {
    std::cerr << argv[0] << ": " << e.what() << '\n';
    Repository::cleanup(); 
    return EXIT_FAILURE;
  }
  catch (...) {
    std::cerr << argv[0] << ": Unknown exception caught.\n";
    Repository::cleanup(); 
    return EXIT_FAILURE;
  }
  
}


void printUsageAndExit() {
  std::cerr << gengetopt_args_info_usage << '\n';
  Repository::cleanup();
  exit( EXIT_FAILURE );
}



void HerwigInit(string infile, string reponame) {
  breakThePEG();
  {
#   ifdef HERWIG_PKGLIBDIR
    DynamicLoader::appendPath(HERWIG_PKGLIBDIR);
#   endif
#   ifdef THEPEG_PKGLIBDIR
    DynamicLoader::appendPath(THEPEG_PKGLIBDIR);
#   endif
    HoldFlag<> setup(InterfaceBase::NoReadOnly);
    Repository::read(infile, cout);
    Repository::update();
  }
  Repository::save(reponame);
}



void HerwigRead(string reponame, string runname) {
#ifdef HERWIG_PKGDATADIR
  ifstream test(reponame.c_str());
  if ( !test ) {
    reponame = string(HERWIG_PKGDATADIR) + '/' + reponame;
  }
  test.close();
#endif
  Repository::load(reponame);
  breakThePEG();
  if ( !runname.empty() && runname != "-" )
    Repository::read(runname, std::cout);
  else
    Repository::read(std::cin, std::cout, "Herwig++> ");
}



void HerwigRun(string runname, int seed, long N, 
	       bool tics, bool resume, bool keepid) {
  PersistentIStream is(runname, keepid);
  ThePEG::EGPtr eg;
  is >> eg;

  // debugging breakpoint
  breakThePEG();

  if ( !eg ) {
    std::cerr << "Herwig++: EventGenerator not available.\n"
	      << "Check if '" << runname << "' is a valid run file.\n";
    Repository::cleanup();
    exit( EXIT_FAILURE );
  }

  if ( seed > 0 ) eg->setSeed(seed);

  eg->go( resume ? -1 : 1, N, tics );

  if ( tics )
    std::cout << '\n';
}

