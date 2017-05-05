// -*- C++ -*-
//
// HerwigCLI.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "HerwigCLI.h"
#include "herwigopts.h"

#include <ThePEG/Utilities/DynamicLoader.h>
#include <ThePEG/Utilities/Debug.h>
#include <ThePEG/Repository/Repository.h>
#include <ThePEG/Handlers/SamplerBase.h>

namespace Herwig {

void HerwigCLI::quitWithHelp() const {
  std::cerr << gengetopt_args_info_usage << '\n';
  quit();
}

void HerwigCLI::quit() const {
  ThePEG::Repository::cleanup();
  exit( EXIT_FAILURE );
}

HerwigCLI::~HerwigCLI() {
  ThePEG::Repository::cleanup();
}

HerwigCLI::HerwigCLI(int argc, char * argv[]) 
  : runMode_(RunMode::ERROR), 
    resume_(false), tics_(true), tag_(),
    inputfile_(), repository_(), setupfile_(),
    integrationList_(),
    N_(-1), seed_(0), jobs_(1),
    jobsize_(0), maxjobs_(0)
{
  gengetopt_args_info args_info;
    
  if ( cmdline_parser( argc, argv, &args_info ) != 0 ) {
    std::cerr << "Could not parse command line.\n";
    return;
  }

  if ( args_info.version_given ) {
    std::cout << 
#include "hgstamp.inc"
"" << '\n';
    std::cout << ThePEG::Repository::version() << std::endl;
    cmdline_parser_free( &args_info );
    exit( EXIT_SUCCESS );
  }

  // require one command
  if ( args_info.inputs_num < 1 )
    quitWithHelp();

  // Define runMode of program
  std::string tmpRunMode = args_info.inputs[0];
  if      ( tmpRunMode == "init" )       { runMode_ = RunMode::INIT; }
  else if ( tmpRunMode == "read" )       { runMode_ = RunMode::READ; }
  else if ( tmpRunMode == "build" )      { runMode_ = RunMode::BUILD; }
  else if ( tmpRunMode == "integrate" )  { runMode_ = RunMode::INTEGRATE; }
  else if ( tmpRunMode == "mergegrids" ) { runMode_ = RunMode::MERGEGRIDS; }
  else if ( tmpRunMode == "run" )        { runMode_ = RunMode::RUN; }
  else {
    runMode_ = RunMode::ERROR;
    quitWithHelp();
  }

  // Use second argument as input- or runfile name
  if ( args_info.inputs_num > 1 )
    inputfile_ = args_info.inputs[1];
  
  // Defaults for these filenames are set in the ggo file
  repository_ = args_info.repo_arg;

  // Number of events
  if ( args_info.numevents_given )
    N_ = args_info.numevents_arg;

  // RNG seed
  if ( args_info.seed_given ) {
    seed_ = args_info.seed_arg;
  }

  // run name tag (default given in ggo file)
  tag_ = args_info.tag_arg;

  // run modification file
  if ( args_info.setupfile_given )
    setupfile_ = args_info.setupfile_arg;

  // parallel jobs
  if ( args_info.jobs_given )
    jobs_ = args_info.jobs_arg;
  
  // Directories from which Herwig reads filesystemfor ( size_t i = 0; i < args_info.append_read_given; ++i )
  for ( size_t i = 0; i < args_info.append_read_given; ++i )
    appendReadDirectories_.push_back( args_info.append_read_arg[i] );
  for ( size_t i = 0; i < args_info.prepend_read_given; ++i )
    prependReadDirectories_.push_back( args_info.prepend_read_arg[i] );

  // Library search path for dlopen()
  for ( size_t i = 0; i < args_info.append_given; ++i )
    ThePEG::DynamicLoader::appendPath( args_info.append_arg[i] );
  for ( size_t i = 0; i < args_info.prepend_given; ++i )
    ThePEG::DynamicLoader::prependPath( args_info.prepend_arg[i] );
  
  // Debugging level
  if ( args_info.debug_given )
    ThePEG::Debug::setDebug( args_info.debug_arg );

  // Floating point exceptions
  if ( args_info.debug_fpe_flag ) 
    ThePEG::Debug::unmaskFpuErrors();

  // Exit-on-error flag
  if ( ! args_info.noexitonerror_flag )
    ThePEG::Repository::exitOnError() = 1;

  // Tics
  if ( args_info.quiet_flag )
    tics_ = false;

  // integration list
  if ( args_info.jobid_given ) {
    integrationList_ = "integrationJob" + std::string(args_info.jobid_arg);
  }

  // job size
  if ( args_info.jobsize_given ) {
    if ( runMode_ != RunMode::BUILD ) {
      std::cerr << "--jobsize option is only available in 'build' mode.\n";
      quitWithHelp();
    }
    jobsize_ = args_info.jobsize_arg;
    ThePEG::SamplerBase::setIntegratePerJob(jobsize_);
  }

  // max integration jobs
  if ( args_info.maxjobs_given ) {
    if ( runMode_ != RunMode::BUILD ) {
      std::cerr << "--maxjobs option is only available in 'build' mode.\n";
      quitWithHelp();
    }
    maxjobs_ = args_info.maxjobs_arg;
    ThePEG::SamplerBase::setIntegrationJobs(maxjobs_);
  }

  // Resume
  if ( args_info.resume_flag )
    resume_ = true;
  
  cmdline_parser_free( &args_info );

}



}
