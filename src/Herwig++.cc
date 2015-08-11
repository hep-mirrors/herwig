// -*- C++ -*-
//
// Herwig++.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "Herwig++.h"

using namespace ThePEG;

int main(int argc, char * argv[]) {
 
  try {

    // read in command line options by using a singleton object
    HelperReadInCommandLineParameters* comLineParam = HelperReadInCommandLineParameters::instance(argc, argv);  
    if(!comLineParam->getReadInWasSuccessful()) {
      std::cerr << "Read in of command line parameters was not successful. Program execution will stop now.";
      return EXIT_FAILURE;
    }
    
  
    // Call program switches according to runMode
    switch ( comLineParam->getRunMode() ) {
    case HerwigRunMode::INIT:        HerwigInit(comLineParam->getRunName(), comLineParam->getRepoName()); break;
    case HerwigRunMode::READ:        HerwigRead(comLineParam->getRepoName(), comLineParam->getRunName()); break;
    case HerwigRunMode::BUILD:       HerwigBuild(comLineParam->getRepoName(), comLineParam->getRunName()); break;
    case HerwigRunMode::INTEGRATE:   HerwigIntegrate(comLineParam->getRunName(), comLineParam->getSetupFile(),
						      comLineParam->getSeed(), comLineParam->getTag(), 
						      comLineParam->getN(), comLineParam->getTics(),
						      comLineParam->getResume(), comLineParam->getJobs(), 
						      comLineParam->getIntegrationList());  
				      break;
    case HerwigRunMode::RUN:         HerwigRun(comLineParam->getRunName(), comLineParam->getSetupFile(),
						      comLineParam->getSeed(), comLineParam->getTag(), 
						      comLineParam->getN(), comLineParam->getTics(),
						      comLineParam->getResume(), comLineParam->getJobs(), 
					              comLineParam->getIntegrationList());  
				      break;
    case HerwigRunMode::ERROR:       std::cerr << "Error during read in of command line parameters. Program execution will stop now."; return EXIT_FAILURE;
    default:          		      printUsageAndExit();
    }


    // Clean repository after program was running.
    Repository::cleanup();     
    return EXIT_SUCCESS;

    // end try block, catching of exceptions
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

void printUsageAndExit() {
  std::cerr << gengetopt_args_info_usage << '\n';
  Repository::cleanup();
  exit( EXIT_FAILURE );
}

void HerwigInit(string infile, string reponame) {
  SamplerBase::setRunLevel(SamplerBase::InitMode);
  
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
  Repository::save(reponame);
}

void HerwigRead(string reponame, string runname) {
  SamplerBase::setRunLevel(SamplerBase::ReadMode);
  
  HerwigGenericRead(reponame, runname);
}

void HerwigBuild(string reponame, string runname) {
  SamplerBase::setRunLevel(SamplerBase::BuildMode);
  
  HerwigGenericRead(reponame, runname);
}

//void HerwigGenericRead(string reponame, string runname, const gengetopt_args_info & args_info)
void HerwigGenericRead(string reponame, string runname) {
#ifdef HERWIG_PKGDATADIR
  ifstream test(reponame.c_str());
  if ( !test ) {
    reponame = string(HERWIG_PKGDATADIR) + '/' + reponame;
  }
  test.close();
#endif
  string msg = Repository::load(reponame);
  if ( ! msg.empty() ) cerr << msg << '\n';
  
  // Invocation of setSearchPaths is necessary since Repository::load revokes all setted paths.  
  // Singleton class is used to return instance
  HelperReadInCommandLineParameters::instance()->setSearchPaths();
  
  breakThePEG();
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
	       int seed, string tag, long int N,
	       bool tics, bool resume, int jobs,
	       string integrationList) {
  SamplerBase::setRunLevel(SamplerBase::RunMode);
  
  HerwigGenericRun(runname, setupfile, seed, tag, N, tics, resume, jobs, false, integrationList);
}

void HerwigIntegrate(string runname, string setupfile,
	       int seed, string tag, long int N,
	       bool tics, bool resume, int jobs,
	       string integrationList) {
  SamplerBase::setRunLevel(SamplerBase::IntegrationMode);
  
  HerwigGenericRun(runname, setupfile, seed, tag, N, tics, resume, jobs, true, integrationList);
}

void HerwigGenericRun(string runname, string setupfile,
	       int seed, string tag, long N, 
	       bool tics, bool resume, int jobs,
	       bool integrationJob,
	       string integrationList) {
  // If runMode is integration or run, we need a runname
  if (runname.empty() ) {
    std::cerr << "Error: You need to supply a runfile name.\n";
    printUsageAndExit();
  }

  if ( integrationJob && jobs > 1 ) {
    std::cerr << "parallel event generation is not applicable to integrate\n";
    exit( EXIT_FAILURE );
  }

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

  Herwig::RunDirectories::pushRunId(eg->runName());
  if ( !setupfile.empty() )
    Herwig::RunDirectories::pushRunId(setupfile);
  if ( !tag.empty() )
    Herwig::RunDirectories::pushRunId(tag);
  if ( !integrationList.empty() )
    Herwig::RunDirectories::pushRunId(integrationList);

  if ( seed > 0 ) {
    ostringstream sseed;
    sseed << seed;
    Herwig::RunDirectories::pushRunId(sseed.str());
  }

  if ( seed > 0 ) eg->setSeed(seed);
  if ( !setupfile.empty() ) eg->addTag("-" + setupfile);
  if ( !tag.empty() ) eg->addTag(tag);

  if ( integrationJob ) {
    Ptr<StandardEventHandler>::tptr eh =
      dynamic_ptr_cast<Ptr<StandardEventHandler>::tptr>(eg->eventHandler());
    if ( !eh ) {
      std::cerr << "Herwig++: Cannot set integration mode for a non-standard EventHandler.\n";
      Repository::cleanup();
      exit( EXIT_FAILURE );
    }
    if ( !integrationList.empty() )
      eh->sampler()->integrationList(integrationList);
  }

  if ( ! setupfile.empty() ) {
    SamplerBase::setupFileUsed();
    string msg = Repository::modifyEventGenerator(*eg, setupfile, cout, true);
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
    pid_t pid;

    const int maxlen = log10(jobs) + 1;

    for (int n=0; n<jobs; n++) {
      ostringstream tmp;
      tmp << std::setfill('0') << std::setw(maxlen) << n+1;
      const string nstr = tmp.str();

      pid = fork();
      if (pid == -1) { // fork failed
        std::perror("Herwig++: fork");
        Repository::cleanup();
        exit( EXIT_FAILURE );
      }
      else if ( pid == 0 ) { // we're the child
        if ( tics ) std::cout << "Forked child " << n << ", PID " << getpid() << std::endl;
        eg->setSeed( seed + n );
        eg->addTag( tag + "-" + nstr );
	Herwig::RunDirectories::pushRunId( nstr );
        eg->go( resume ? -1 : 1, N / jobs, false );
        break; // avoid sub-forks
      }
      // nothing to do here if we're the parent
    }

    // children have nothing else to do
    if (pid == 0) return;

    if ( tics ) std::cout << "Waiting for forked jobs." << std::endl;
    int status;
    pid_t child;
    bool cleanrun = true;
    while (true) {
    	child = wait(&status);
    	if (child == -1) {
    		if (errno == ECHILD) {
    			if ( tics ) std::cout << "No more forked jobs." << std::endl;
    			break;
    		}
    		else {
	                std::perror("Herwig++: waitpid");
	                Repository::cleanup();
        	        exit(EXIT_FAILURE);
        	}
        }

        if (WIFEXITED(status)) {
            if (WEXITSTATUS(status) != 0) 
            	cleanrun = false;
            if ( tics ) std::cout << "PID " << child 
                      << " exited, status=" << WEXITSTATUS(status)
                      << std::endl;
        } else if (WIFSIGNALED(status)) {
            // a clean SIGTERM is handled in the child
            // and will count as exit above, so...
            cleanrun = false;
            if ( tics ) std::cout << "PID " << child 
           	<< " killed by signal " << WTERMSIG(status)
           	<< std::endl;
        }
    }
    if (! cleanrun) {
    	Repository::cleanup();
    	exit(EXIT_FAILURE);
    }
  }
}

HelperReadInCommandLineParameters::HelperReadInCommandLineParameters(int argc, char* argv[]) {
  // Define default values for errors first
  m_readInWasSuccessful = false;
  m_runMode = HerwigRunMode::ERROR;
  
  
  // read command line options 
    
  if ( cmdline_parser( argc, argv, &m_args_info ) != 0 ) {
    std::cerr << "Could not parse command line.\n";
    return;
  }

  // require one command
  if ( m_args_info.inputs_num < 1 )
    printUsageAndExit();

  // Define runMode of program
  std::string tmpRunMode = m_args_info.inputs[0];
  if ( tmpRunMode == "init" ) {
    m_runMode = HerwigRunMode::INIT;
  } else if ( tmpRunMode == "read" ) {
    m_runMode = HerwigRunMode::READ;
  } else if ( tmpRunMode == "build" ) {
    m_runMode = HerwigRunMode::BUILD;
  } else if ( tmpRunMode == "integrate" ) {
    m_runMode = HerwigRunMode::INTEGRATE;
  } else if ( tmpRunMode == "run" ) {
    m_runMode = HerwigRunMode::RUN;
  } else {
    m_runMode = HerwigRunMode::ERROR;
    printUsageAndExit();
  }

  // Use second argument as input- or runfile name
  if ( m_args_info.inputs_num > 1 )
    m_runname = m_args_info.inputs[1];
  
  // Defaults for these filenames are set in the ggo file
  m_reponame = m_args_info.repo_arg;

  // Number of events
  m_N = -1;
  if ( m_args_info.numevents_given )
    m_N = m_args_info.numevents_arg;

  // RNG seed
  m_seed = 0;
  if ( m_args_info.seed_given ) {
    m_seed = m_args_info.seed_arg;
  }

  // run name tag (default given in ggo file)
  m_tag = m_args_info.tag_arg;

  // run modification file
  m_setupfile = "";
  if ( m_args_info.setupfile_given )
    m_setupfile = m_args_info.setupfile_arg;

  // parallel jobs
  m_jobs = 1;
  if ( m_args_info.jobs_given )
    m_jobs = m_args_info.jobs_arg;

  this->setSearchPaths();

  // Library search path for dlopen()
  for ( size_t i = 0; i < m_args_info.append_given; ++i )
    DynamicLoader::appendPath( m_args_info.append_arg[i] );
  for ( size_t i = 0; i < m_args_info.prepend_given; ++i )
    DynamicLoader::prependPath( m_args_info.prepend_arg[i] );
  
  // Debugging level
  if ( m_args_info.debug_given )
    Debug::setDebug( m_args_info.debug_arg );

  // Floating point exceptions
  if ( m_args_info.debug_fpe_flag ) 
    Debug::unmaskFpuErrors();

  // Exit-on-error flag
  if ( ! m_args_info.noexitonerror_flag )
    Repository::exitOnError() = 1;

  // Tics
  m_tics = true;
  if ( m_args_info.quiet_flag )
    m_tics = false;

  // integration list
  m_integrationList = "";
  if ( m_args_info.jobid_given ) {
    m_integrationList = "integrationJob" + string(m_args_info.jobid_arg);
  }

  // job size
  m_jobsize = 0;
  if ( m_args_info.jobsize_given ) {
    m_jobsize = m_args_info.jobsize_arg;
    SamplerBase::setIntegratePerJob(m_jobsize);
  }

  // max integration jobs
  m_maxjobs = 0;
  if ( m_args_info.maxjobs_given ) {
    m_maxjobs = m_args_info.maxjobs_arg;
    SamplerBase::setIntegrationJobs(m_maxjobs);
  }

  // Resume
  m_resume = false;
  if ( m_args_info.resume_flag )
    m_resume = true;
  
  m_readInWasSuccessful = true;
  return;
}

void HelperReadInCommandLineParameters::setSearchPaths() {
  // Search path for read command uses CWD first
  string cwd = boost::filesystem::current_path().string();
  Repository::prependReadDir( cwd );
  // append command line choices
  for ( size_t i = 0; i < m_args_info.append_read_given; ++i )
    Repository::appendReadDir( m_args_info.append_read_arg[i] );
  for ( size_t i = 0; i < m_args_info.prepend_read_given; ++i )
    Repository::prependReadDir( m_args_info.prepend_read_arg[i] );
}

