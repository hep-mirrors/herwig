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
#include <ThePEG/Handlers/StandardEventHandler.h>
#include <ThePEG/Handlers/SamplerBase.h>
#include <iostream>

#include <config.h>
#ifdef HAVE_UNISTD_H
#include <queue>
#include <unistd.h>
#include <sys/wait.h>
#endif


using namespace ThePEG;

void printUsageAndExit();

void HerwigInit(string infile, string reponame);
void HerwigRead(string reponame, string runname,
		const gengetopt_args_info & args_info,
		bool postponeInitialize);
void HerwigRun(string runname, string setupfile,
	       int seed, string tag, long N, 
	       bool tics, bool resume, int jobs,
	       bool integrationJob,
	       string integrationList);

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
    enum { INIT, READ, BUILD, INTEGRATE, RUN, ERROR } status;
    std::string runType = args_info.inputs[0];
    if ( runType == "init" )
      status = INIT;
    else if ( runType == "read" ) 
      status = READ;
    else if ( runType == "build" ) 
      status = BUILD;
    else if ( runType == "integrate" ) 
      status = INTEGRATE;
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
    if ( ( status == RUN || status == INTEGRATE ) && runname.empty() ) {
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

    // run modifccation file
    string setupfile = "";
    if ( args_info.setupfile_given )
      setupfile = args_info.setupfile_arg;

    // parallel jobs
    int jobs = 1;
#   ifdef HAVE_UNISTD_H
    if ( args_info.jobs_given )
      jobs = min( args_info.jobs_arg, 10 );
#   endif

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

    // integration list
    string integrationList = "";
    if ( args_info.bins_given )
      integrationList = args_info.bins_arg;

    // Resume
    bool resume = false;
    if ( args_info.resume_flag )
      resume = true;

    // *** End of command line parsing ***
   
    // Call mode
    switch ( status ) {
    case INIT:        HerwigInit( runname, reponame ); break;
    case READ:        HerwigRead( reponame, runname, args_info, false ); break;
    case BUILD:       HerwigRead( reponame, runname, args_info, true ); break;
    case INTEGRATE:   HerwigRun( runname, setupfile , seed, tag, N, tics, resume, jobs, true, integrationList );  break;
    case RUN:         HerwigRun( runname, setupfile , seed, tag, N, tics, resume, jobs, false, integrationList );  break;
    default:          printUsageAndExit();
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
		const gengetopt_args_info & args_info,
		bool postponeInitialize) {
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
  if ( postponeInitialize )
    SamplerBase::doPostponeInitialize();
  if ( !runname.empty() && runname != "-" ) {
    string msg = Repository::read(runname, std::cout);
    if ( ! msg.empty() ) cerr << msg << '\n';
  }
  else {
    Repository::exitOnError() = 0;
    Repository::read(std::cin, std::cout, "Herwig++> ");
  }
}



void HerwigRun(string runname, string setupfile,
	       int seed, string tag, long N, 
	       bool tics, bool resume, int jobs,
	       bool integrationJob,
	       string integrationList) {
  PersistentIStream is(runname);
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

  if ( integrationJob ) {
    Ptr<StandardEventHandler>::tptr eh =
      dynamic_ptr_cast<Ptr<StandardEventHandler>::tptr>(eg->eventHandler());
    if ( !eh ) {
      std::cerr << "Herwig++: Cannot set integration mode for a non-standard EventHandler.\n";
      Repository::cleanup();
      exit( EXIT_FAILURE );
    }
    eh->sampler()->isIntegrationJob();
    if ( integrationList != "" )
      eh->sampler()->integrationList(integrationList);
  }

  if ( ! setupfile.empty() ) {
    string msg = Repository::modifyEventGenerator(*eg, setupfile, cout, integrationJob);
    if ( ! msg.empty() ) cerr << msg << '\n';
    if ( integrationJob )
      return;
  }

  if ( integrationJob ) {
    Repository::resetEventGenerator(*eg);
    return;
  }
  
  if (jobs <= 1) {

    eg->go( resume ? -1 : 1, N, tics );
    if ( tics ) std::cout << '\n';
  
  }
  else { // forked jobs

#   ifdef HAVE_UNISTD_H

    std::queue<pid_t> pids;
    pid_t pid;

    for (int n=0; n<jobs; n++) {
      pid = fork();
      if (pid == -1) {
        std::cerr << "Herwig++: Problem in fork().\n";
        Repository::cleanup();
        exit( EXIT_FAILURE );
      }
      else if ( pid == 0 ) {
        std::cout << "Forked child " << n << ", PID " << getpid() << std::endl;
        eg->setSeed( seed + n );
        // fix numbering to allow n > 10
        assert( n <= 10 );
        eg->addTag( tag + "-" + char( 48 + n ) );
        eg->go( resume ? -1 : 1, N / jobs, false );
        break;
      }
      else {
        pids.push(pid);
      }
    }

    if (pid == 0) return;

    while (! pids.empty() ) {
      std::cout << "Waiting for " << pids.size() << " job(s)." << std::endl;
      waitpid(pids.front(), NULL, 0);
      std::cout << "PID " << pids.front() << " done." << std::endl;
      pids.pop();
    }

#   endif
    return;

  }
}
