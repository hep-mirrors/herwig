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
void HerwigRead(string reponame, string runname,
		const gengetopt_args_info & args_info);
void HerwigRun(string runname, int seed, string tag, long N, 
	       bool tics, bool resume, bool keepid);

void setSearchPaths(const gengetopt_args_info & args_info);

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
    enum { INIT, READ, RUN, ERROR } status;
    std::string runType = args_info.inputs[0];
    if ( runType == "init" )
      status = INIT;
    else if ( runType == "read" ) 
      status = READ;
    else if ( runType == "run" )  
      status = RUN;
    else {
      status = ERROR;
      printUsageAndExit();
    }

    // Use second argument as input- or runfile name
    string runname;
    if ( args_info.inputs_num > 1 )
      runname = args_info.inputs[1];

    // If status is RUN, we need a runname
    if ( status == RUN && runname.empty() ) {
      cerr << "Error: You need to supply a runfile name.\n";
      printUsageAndExit();
    }

    // If status is INIT, we need a runname
    if ( status == INIT && runname.empty() )
	runname = "HerwigDefaults.in";

    // Defaults for these filenames are set in the ggo file
    std::string reponame = args_info.repo_arg;

    // Number of events
    long N = -1;
    if ( args_info.numevents_given )
      N = args_info.numevents_arg;

    // RNG seed
    int seed = 0;
    if ( args_info.seed_given )
      seed = args_info.seed_arg;

    // run name tag (default given in ggo file)
    string tag = args_info.tag_arg;

    setSearchPaths(args_info);
  
    // Library search path for dlopen()
    for ( size_t i = 0; i < args_info.append_given; ++i )
      DynamicLoader::appendPath( args_info.append_arg[i] );
    for ( size_t i = 0; i < args_info.prepend_given; ++i )
      DynamicLoader::prependPath( args_info.prepend_arg[i] );
    
    // Debugging level
    if ( args_info.debug_given )
      Debug::setDebug( args_info.debug_arg );

    // Floating point exceptions
    if ( args_info.debug_fpe_flag ) 
      Debug::unmaskFpuErrors();

    // Exit-on-error flag
    if ( ! args_info.noexitonerror_flag )
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
    case INIT:  HerwigInit( runname, reponame ); break;
    case READ:  HerwigRead( reponame, runname, args_info ); break;
    case RUN:   HerwigRun( runname, seed, tag, N, 
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
  catch (const char* what) {
    std::cerr << argv[0] << ": caught exception: "
	      << what << "\n";
    Repository::cleanup(); 
    return EXIT_FAILURE;
  }
  catch (...) {
    std::cerr << argv[0] << ": Unknown exception caught.\n";
    Repository::cleanup(); 
    return EXIT_FAILURE;
  }
  
}


void setSearchPaths(const gengetopt_args_info & args_info) {
  // Search path for read command
  for ( size_t i = 0; i < args_info.append_read_given; ++i )
    Repository::appendReadDir( args_info.append_read_arg[i] );
  for ( size_t i = 0; i < args_info.prepend_read_given; ++i )
    Repository::prependReadDir( args_info.prepend_read_arg[i] );
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
    string msg = Repository::read(infile, cout);
    if ( ! msg.empty() ) cerr << msg << '\n';
    Repository::update();
  }
  Repository::save(reponame);
}


void HerwigRead(string reponame, string runname, 
		const gengetopt_args_info & args_info) {
#ifdef HERWIG_PKGDATADIR
  ifstream test(reponame.c_str());
  if ( !test ) {
    reponame = string(HERWIG_PKGDATADIR) + '/' + reponame;
  }
  test.close();
#endif
  string msg = Repository::load(reponame);
  if ( ! msg.empty() ) cerr << msg << '\n';
  setSearchPaths(args_info);
  breakThePEG();
  if ( !runname.empty() && runname != "-" ) {
    string msg = Repository::read(runname, std::cout);
    if ( ! msg.empty() ) cerr << msg << '\n';
  }
  else
    Repository::read(std::cin, std::cout, "Herwig++> ");
}



void HerwigRun(string runname, int seed, string tag, long N, 
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
  if ( !tag.empty() ) eg->addTag(tag);

  eg->go( resume ? -1 : 1, N, tics );

  if ( tics )
    std::cout << '\n';
}

