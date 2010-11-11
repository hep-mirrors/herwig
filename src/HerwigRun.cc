// -*- C++ -*-
//
// HerwigRun.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "HerwigRun.h"
#include "versionstring.h"
#include <ThePEG/Utilities/DynamicLoader.h>
#include <ThePEG/Utilities/Debug.h>
#include <ThePEG/Repository/Repository.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/EventRecord/Event.h>
#include <cassert>
#include <iostream>

using namespace Herwig;
using namespace ThePEG;

const string usage = " init|read|run "
                   " [-N num-events] [--seed random-generator-seed] "
                   "[-d debug-level] "
                   "[--exitonerror] [-l load-path] "
                   "[-L first-load-path] [-r repo-file] "
                   "[-i initialization file] [run-file]\n";

HerwigRun::HerwigRun(int argc, char **argv) 
: N_(-1),
  seed_(0),
  runname_(),
  status_(UNKNOWN)
{
  std::string reponame("HerwigDefaults.rpo");
  std::string infile("HerwigDefaults.in");

  std::string runType;
  if ( argc > 1 ) 
    runType = argv[1];

  if ( runType == "init" ) 
    status_ = INIT;
  else if ( runType == "read" ) 
    status_ = READ;
  else if ( runType == "run" ) 
    status_ = RUN;
  else if ( runType == "-h" || runType == "--help" ) {
    std::cout << "Usage: " << argv[0] << usage;
    printHelp(std::cout);
    status_ = SUCCESS;
    return;
  } 
  else if ( runType == "-v" || runType == "--version" ) {
    std::cout << versionstring << std::endl;
    status_ = SUCCESS;
    return;
  }
  else {
    std::cerr << "Usage: " << argv[0] << usage;
    status_ = ERROR;
    return;
  }


  for ( int iarg = 2; iarg < argc; ++iarg ) {
    std::string arg = argv[iarg];
    if ( arg == "-r" ) 
      reponame = argv[++iarg];
    else if( arg == "-i" ) 
      infile = argv[++iarg];
    else if ( arg == "-N" ) 
      N_ = atoi(argv[++iarg]);
    else if ( arg.substr(0,2) == "-N" ) 
      N_ = atoi(arg.substr(2).c_str());
    else if ( arg == "-seed" || arg == "--seed" ) 
      seed_ = atoi(argv[++iarg]);
    else if ( arg == "-l" ) 
      DynamicLoader::appendPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-l" )
      DynamicLoader::appendPath(arg.substr(2));
    else if ( arg == "-L" ) 
      DynamicLoader::prependPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-L" )
      DynamicLoader::prependPath(arg.substr(2));
    else if ( arg == "-d" ) 
      Debug::setDebug(atoi(argv[++iarg]));
    else if ( arg.substr(0,2) == "-d" )
      Debug::setDebug(atoi(arg.substr(2).c_str()));
    else if ( arg == "--exitonerror" ) 
      Repository::exitOnError() = 1;
    else if ( arg == "-h" || arg == "--help" ) {
      std::cout << "Usage: " << argv[0] << usage;
      printHelp(std::cout);
      status_ = SUCCESS;
      return;
    }
    else
      runname_ = arg;
  }
  if ( Debug::level ) 
    Debug::unmaskFpuErrors();


  if ( status_ == INIT ) {
    // debugging breakpoint
    breakThePEG();
    {
#ifdef HERWIG_PKGLIBDIR
      DynamicLoader::appendPath(HERWIG_PKGLIBDIR);
#endif
#ifdef THEPEG_PKGLIBDIR
      DynamicLoader::appendPath(THEPEG_PKGLIBDIR);
#endif
      HoldFlag<> setup(InterfaceBase::NoReadOnly);
      Repository::read(infile, cout);
      Repository::update();
    }
    Repository::save(reponame);
    status_ = SUCCESS;
  } 

  else if ( status_ == READ ) {
#ifdef HERWIG_PKGDATADIR
    ifstream test(reponame.c_str());
    if ( !test ) {
      reponame = string(HERWIG_PKGDATADIR) + '/' + reponame;
    }
    test.close();
#endif
    Repository::load(reponame);
    breakThePEG();
    if ( !runname_.empty() && runname_ != "-" ) {
      Repository::read(runname_, std::cout);
    } else {
      Repository::read(std::cin, std::cout, "Herwig++> ");
    }
    status_ = SUCCESS;
  } 
  else if ( runname_.empty() ) {
    std::cerr << "No run-file specified.\n";
    status_ = ERROR;
  } 
  else if ( status_ == RUN ) {
    generateEvents();
  }
  else {
    std::cerr << "Argument parse error.\n"
	      << "Usage: " << argv[0] << usage;
    status_ = ERROR;
  }
}
  

void HerwigRun::generateEvents() {
  assert( status_ == RUN );

  PersistentIStream is(runname_);
  ThePEG::EGPtr eg;
  is >> eg;

  // debugging breakpoint
  breakThePEG();

  if ( !eg ) {
    std::cerr << __FILE__ << ": EventGenerator not available.\n";
    status_ = ERROR;
    return;
  }
  if ( seed_ > 0 ) 
    eg->setSeed(seed_);
  eg->go( 1, N_, true );
  std::cout << '\n';
  status_ = SUCCESS;
}

bool HerwigRun::good() const { return status_ == SUCCESS; }

void HerwigRun::printHelp(std::ostream &out) {
  out << '\n'
      << "One of the following options is required.\n"
      << "==============================================================\n"
      << "init    - Reread default file and create .rpo\n"
      << "read    - Read input file and create run file\n"
      << "run     - Read run file and run\n"
      << "==============================================================\n"
      << "These are optional commands\n"
      << "==============================================================\n"
      << "-N      - Set number of events in the run\n"
      << "-d      - Sets the ThePEG debug level (see ThePEG::Debug)\n"
      << "-r      - Changes the repository file from HerwigDefaults.rpo\n"
      << "-i      - Changes the repo input file from HerwigDefaults.in\n"
      << "-l      - Adds path to dynamically load library (to end)\n"
      << "-L      - Adds path to dynamically load library (to beginning)\n"
      << "-seed   - Sets the random seed on initialization\n"
      << "-h      - Displays this help message\n"
      << "==============================================================\n";
}

HerwigRun::~HerwigRun() { 
  Repository::cleanup(); 
}
