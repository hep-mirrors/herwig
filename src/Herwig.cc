// -*- C++ -*-
//
// Herwig.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "Herwig.h"
#include "HerwigUI.h"

#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/EventGenerator.h>

#include <ThePEG/Utilities/DynamicLoader.h>
#include <ThePEG/Utilities/Debug.h>
#include <ThePEG/Repository/Repository.h>
#include <ThePEG/Handlers/SamplerBase.h>

#include <ThePEG/Handlers/StandardEventHandler.h>

#include <iostream>
#include <sstream>

#include <unistd.h>
#include <sys/wait.h>

#include <boost/filesystem.hpp>

#include "Herwig/Utilities/RunDirectories.h"

using namespace ThePEG;

namespace {
/**
 * Search paths for Repository read
 * 
 * You can pass two string vectors with directories which Herwig will use to look in for files.
 * A vector with directories which will be prepended and a vector with directories which will be appended.
 * Both vectors are optional.
 */
void setSearchPaths(const Herwig::HerwigUI & ui,
                    bool usePWD = true);

void HerwigGenericRead(const Herwig::HerwigUI & ui);
void HerwigGenericRun(const Herwig::HerwigUI & ui);
}

namespace Herwig {
namespace API {

void init(const HerwigUI & ui) {
  setSearchPaths(ui, false);
  SamplerBase::setRunLevel(SamplerBase::InitMode);
  
  string infile = ui.inputfile();

  // If status is INIT, we need an infile name / runname
  if (infile.empty())
    infile = "HerwigDefaults.in";
  
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
  Repository::save(ui.repository());
}

void read(const HerwigUI & ui) {
  setSearchPaths(ui);
  SamplerBase::setRunLevel(SamplerBase::ReadMode);
  HerwigGenericRead(ui);
}

void build(const HerwigUI & ui) {
  setSearchPaths(ui);
  SamplerBase::setRunLevel(SamplerBase::BuildMode);
  HerwigGenericRead(ui);
}

void integrate(const HerwigUI & ui) {
  setSearchPaths(ui);
  SamplerBase::setRunLevel(SamplerBase::IntegrationMode);
  HerwigGenericRun(ui);
}

void run(const HerwigUI & ui) {
  setSearchPaths(ui);
  SamplerBase::setRunLevel(SamplerBase::RunMode);  
  HerwigGenericRun(ui);
}

}
}


namespace {
void setSearchPaths(const Herwig::HerwigUI & ui,
                    bool usePWD) {
  // Search path for read command uses CWD first
  if ( usePWD ) {
        string cwd = boost::filesystem::current_path().string();
        Repository::prependReadDir( cwd );
  }
  // append command line choices for directories from which Herwig will read.
  Repository::appendReadDir(ui.appendReadDirectories());
  Repository::prependReadDir(ui.prependReadDirectories());
}

//void HerwigGenericRead(string reponame, string runname, const gengetopt_args_info & args_info)
void HerwigGenericRead(const Herwig::HerwigUI & ui) {
  string reponame = ui.repository();
#ifdef HERWIG_PKGDATADIR
  ifstream test(reponame.c_str());
  if ( !test ) {
    reponame = string(HERWIG_PKGDATADIR) + '/' + reponame;
  }
  test.close();
#endif
  string msg = Repository::load(reponame);
  if ( ! msg.empty() ) cerr << msg << '\n';
  
  // Repeated invocation of setSearchPaths is necessary since Repository::load revokes all setted paths.  
  setSearchPaths(ui);
  
  breakThePEG();
  if ( !ui.inputfile().empty() && ui.inputfile() != "-" ) {
    string msg = Repository::read(ui.inputfile(), std::cout);
    if ( ! msg.empty() ) cerr << msg << '\n';
  }
  else {
    Repository::exitOnError() = 0;
    Repository::read(std::cin, std::cout, "Herwig> ");
  }
}

void HerwigGenericRun(const Herwig::HerwigUI & ui) {
  // If runMode is integration or run, we need a runname
  const string runname = ui.inputfile();
  if (runname.empty() ) {
    std::cerr << "Error: You need to supply a runfile name.\n";
    ui.quitWithHelp();
  }

  if ( ui.integrationJob() && ui.jobs() > 1 ) {
    std::cerr << "parallel event generation is not applicable to integrate\n";
    ui.quit();
  }

  PersistentIStream is(runname);
  ThePEG::EGPtr eg;
  is >> eg;

  // debugging breakpoint
  breakThePEG();

  if ( !eg ) {
    std::cerr << "Herwig: EventGenerator not available.\n"
	      << "Check if '" << runname << "' is a valid run file.\n";
    ui.quit();
  }

  Herwig::RunDirectories::pushRunId(eg->runName());
  if ( !ui.setupfile().empty() )
    Herwig::RunDirectories::pushRunId(ui.setupfile());
  if ( !ui.tag().empty() )
    Herwig::RunDirectories::pushRunId(ui.tag());
  if ( !ui.integrationList().empty() )
    Herwig::RunDirectories::pushRunId(ui.integrationList());

  if ( ui.seed() > 0 ) {
    ostringstream sseed;
    sseed << ui.seed();
    Herwig::RunDirectories::pushRunId(sseed.str());
  }

  if ( ui.seed() > 0 ) eg->setSeed(ui.seed());
  if ( !ui.setupfile().empty() ) eg->addTag("-" + ui.setupfile());
  if ( !ui.tag().empty() ) eg->addTag("-" + ui.tag());

  if ( ui.integrationJob() ) {
    Ptr<StandardEventHandler>::tptr eh =
      dynamic_ptr_cast<Ptr<StandardEventHandler>::tptr>(eg->eventHandler());
    if ( !eh ) {
      std::cerr << "Herwig: Cannot set integration mode for a non-standard EventHandler.\n";
      ui.quit();
    }
    if ( !ui.integrationList().empty() )
      eh->sampler()->integrationList(ui.integrationList());
  }

  if ( ! ui.setupfile().empty() ) {
    SamplerBase::setupFileUsed();
    string msg = Repository::modifyEventGenerator(*eg, ui.setupfile(), cout, true);
    if ( ! msg.empty() ) cerr << msg << '\n';
    if ( ui.integrationJob() )
      return;
  }

  if ( ui.integrationJob() ) {
    Repository::resetEventGenerator(*eg);
    return;
  }
  
  if (ui.jobs() <= 1) {

    eg->go( ui.resume() ? -1 : 1, ui.N(), ui.tics() );
    if ( ui.tics() ) std::cout << '\n';
  
  }
  else { // forked jobs
    pid_t pid = 0;

    const int maxlen = log10(ui.jobs()) + 1;

    // make sure the parent got initialized once so e.g. grid combination is
    // working
    eg->initialize(true);

    for (int n=0; n<ui.jobs(); n++) {
      ostringstream tmp;
      tmp << std::setfill('0') << std::setw(maxlen) << n+1;
      const string nstr = tmp.str();

      pid = fork();
      if (pid == -1) { // fork failed
        std::perror("Herwig: fork");
        ui.quit();
      }
      else if ( pid == 0 ) { // we're the child
        if ( ui.tics() ) std::cout << "Forked child " << n << ", PID " << getpid() << std::endl;
        eg->setSeed( ui.seed() + n );
        eg->addTag( "-" + nstr );
	Herwig::RunDirectories::pushRunId( nstr );
        eg->go( ui.resume() ? -1 : 1, ui.N() / ui.jobs(), false );
        break; // avoid sub-forks
      }
      // nothing to do here if we're the parent
    }

    // children have nothing else to do
    if (pid == 0) return;

    if ( ui.tics() ) std::cout << "Waiting for forked jobs." << std::endl;
    int status;
    pid_t child;
    bool cleanrun = true;
    while (true) {
    	child = wait(&status);
    	if (child == -1) {
    		if (errno == ECHILD) {
    			if ( ui.tics() ) std::cout << "No more forked jobs." << std::endl;
    			break;
    		}
    		else {
	                std::perror("Herwig: waitpid");
	                ui.quit();
        	}
        }

        if (WIFEXITED(status)) {
            if (WEXITSTATUS(status) != 0) 
            	cleanrun = false;
            if ( ui.tics() ) std::cout << "PID " << child 
                      << " exited, status=" << WEXITSTATUS(status)
                      << std::endl;
        } else if (WIFSIGNALED(status)) {
            // a clean SIGTERM is handled in the child
            // and will count as exit above, so...
            cleanrun = false;
            if ( ui.tics() ) std::cout << "PID " << child 
           	<< " killed by signal " << WTERMSIG(status)
           	<< std::endl;
        }
    }
    if (! cleanrun) {
    	ui.quit();
    }
  }
}
}


