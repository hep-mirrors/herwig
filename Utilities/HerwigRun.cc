// -*- C++ -*-
//
// HerwigRun.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "HerwigRun.h"
#include "HerwigVersion.h"
#include "HwDebug.h"
#include <ThePEG/Utilities/DynamicLoader.h>
#include <ThePEG/Utilities/SystemUtils.h>
#include <ThePEG/Utilities/StringUtils.h>
#include <ThePEG/Repository/Repository.h>
#include <ThePEG/Persistency/PersistentIStream.h>

using namespace Herwig;
using namespace ThePEG;

const string usage = " init|read|run "
                   " [-N num-events] [--seed random-generator-seed] "
                   "[-d debug-level] [-dHw herwig-debug-level] "
                   "[--exitonerror] [-l load-path] "
                   "[-L first-load-path] [-r repo-file] "
                   "[-i initialization file] [run-file]\n";

HerwigRun::HerwigRun(int argc, char **argv) 
: N(-1),
  ngen(0),
  seed(0),
  repo(),
  repoin("HerwigDefaults.in"),
  egCreated(false), 
  Status(UNKNOWN), 
  isInitialized(false)
{
  std::string runType;
  if(argc > 1) runType = argv[1];
  else runType = "";

  if(runType == "init") Status = INIT;
  else if(runType == "read") Status = READ;
  else if(runType == "run") Status = RUN;
  else if(runType == "-h" || runType == "--help") {
    std::cout << "Usage: " << argv[0] << usage;
    printHelp(std::cout);
    Status = SUCCESS;
    return;
  } else {
    std::cerr << "Usage: " << argv[0] << usage;
    Status = ERROR;
    return;
  }

  for ( int iarg = 2; iarg < argc; ++iarg ) {
    std::string arg = argv[iarg];
    if ( arg == "-r" ) repo = argv[++iarg];
    else if(arg == "-i") repoin = argv[++iarg];
    else if ( arg == "-l" ) DynamicLoader::appendPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-l" )
      DynamicLoader::appendPath(arg.substr(2));
    else if ( arg == "-L" ) DynamicLoader::prependPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-L" )
      DynamicLoader::prependPath(arg.substr(2));
    else if ( arg == "-d" ) Debug::level = atoi(argv[++iarg]);
    else if ( arg.substr(0,2) == "-d" && arg.substr(0,4) != "-dHw" )
	Debug::level = atoi(arg.substr(2).c_str());
    else if ( arg == "-dHw" ) Herwig::HwDebug::level = atoi(argv[++iarg]);
    else if ( arg.substr(0,4) == "-dHw" )
	Herwig::HwDebug::level = atoi(arg.substr(4).c_str());
    else if ( arg == "-N" ) N = atoi(argv[++iarg]);
    else if ( arg.substr(0,2) == "-N" ) N = atoi(arg.substr(2).c_str());
    else if ( arg == "-seed" || arg == "--seed" ) seed = atoi(argv[++iarg]);
    else if ( arg == "--exitonerror" ) Repository::exitOnError() = 1;
    else if ( arg == "-h" || arg == "--help" ) {
      std::cout << "Usage: " << argv[0] << usage;
      printHelp(std::cout);
      Status = SUCCESS;
      return;
    }
    else
      run = arg;
  }
  if ( Status == INIT ) {
    // debugging breakpoint
    breakThePEG();
    {
      HoldFlag<> setup(InterfaceBase::NoReadOnly);
      Repository::read(repoin, cout);
      Repository::update();
    }
    if ( repo.empty() )
      repo = "HerwigDefaults.rpo";
    Repository::save(repo);
    Status = SUCCESS;
  } 
  else if ( Status == READ ) {
    if ( repo.empty() ) {
      repo = "HerwigDefaults.rpo";
      ifstream test(repo.c_str());
      if (!test) {
	repo = HerwigVersion::pkgdatadir + "/HerwigDefaults.rpo";
      }
      test.close();
    }
    Repository::load(repo);
    breakThePEG();
    if ( run.size() && run != "-" ) {
      ifstream is(run.c_str());
      Repository::read(is, std::cout);
    } else {
      Repository::read(std::cin, std::cout, "Herwig++> ");
    }
    Status = SUCCESS;
  } 
  else if ( run.empty() ) {
    std::cerr << "No run-file specified.\n";
    Status = ERROR;
  } 
  else if ( Status == RUN ) {
    generateEvents();
  }
  else {
    std::cerr << "Argument parse error.\n";
    Status = ERROR;
  }
}
  
EGPtr HerwigRun::eventGenerator() {
  if( Status != RUN ) 
    return EGPtr();
  if(!egCreated) {
    PersistentIStream is(run);
    is >> eg;
    breakThePEG();
    egCreated = true;
    if ( eg ) 
      {
	if ( seed > 0 ){
	  eg->setSeed(seed);
	}
	if(!isInitialized) {
	  eg->initialize();
	  isInitialized = true;
	}
      }
    return eg;
  } else 
    {
      return eg; 
    }
}

EventPtr HerwigRun::generateEvent() {
  lastEvent = EventPtr();
  if( Status != RUN ) 
    return EventPtr();
  if(!isInitialized) {
    eg->initialize();
    isInitialized = true;
  }
  if(!egCreated) eventGenerator();
  if(eg) {
    if(ngen < N && ngen < eg->N()) { 
      ngen++; 
      try {
	lastEvent = eg->shoot(); 
      }
      catch ( ... ) {
	eg->finish();
	throw;
      }
      return lastEvent; 
    }
    else {
      return EventPtr();
    }
  } else return lastEvent;
}

void HerwigRun::generateEvents() {
  if ( !isRunMode() || !preparedToRun() ) {
    std::cerr << "Error: EventGenerator not available.\n";
    Status = ERROR;
    return;
  }

  if ( getN() > eventGenerator()->N() )
    std::cerr << "Warning: will only generate " 
	      << eventGenerator()->N() << " events;\n"
	      << "Warning: you can increase NumberOfEvents "
	      << "in the input files.\n";
  long number = std::min( getN(), eventGenerator()->N() );
  long step = std::max( number/100l, 1l );
  std::cout << "Generating events.\r" << std::flush;
  for( long i = 1; i <= number; ++i ) {
    generateEvent();
    if ( i % step == 0 )
      std::cout << "Generated event: " << i 
		<< " of " << number << "\r" << std::flush;
  }
  std::cout << '\n';
  eventGenerator()->finalize();
  Status = SUCCESS;
}

long HerwigRun::getN() const { return N; }
long HerwigRun::getNGen() const { return ngen; }
HerwigRun::RunStatus HerwigRun::status() const { return Status; }
bool HerwigRun::good() const { return status() == SUCCESS; }
bool HerwigRun::isRunMode() const { return status() == RUN; }
bool HerwigRun::isReadMode() const { return status() == READ; }
bool HerwigRun::isInitMode() const { return status() == INIT; }

void HerwigRun::printHelp(std::ostream &out) {
  out << endl
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
      << "-dHw    - Sets the Herwig debug level (see Herwig::Debug)\n"
      << "-r      - Changes the repository file from HerwigDefaults.rpo\n"
      << "-i      - Changes the repo input file from HerwigDefaults.in\n"
      << "-l      - Adds path to dynamically load library (to end)\n"
      << "-L      - Adds path to dynamically load library (to beginning)\n"
      << "-seed   - Sets the random seed on initialization\n"
      << "-h      - Displays this help message\n"
      << "==============================================================\n";
}

bool HerwigRun::preparedToRun() {
  if(eventGenerator()) return true;
  else return false;
}
