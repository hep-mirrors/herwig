// -*- C++ -*-
//
// Herwig++.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef SRC_HERWIG_H
#define SRC_HERWIG_H

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
#include <sstream>
#include <cstdio>
#include <unistd.h>
#include <sys/wait.h>
#include <boost/filesystem.hpp>

#include "Herwig++/Utilities/RunDirectories.h"

using namespace ThePEG;

/**
 * Show usage information in terminal
 * 
 */
void printUsageAndExit();

/**
 * Herwig init run mode
 * 
 */
void HerwigInit(string infile, string reponame);

/**
 * Herwig generic read mode
 * 
 * Function is used by built and read mode of Herwig++.
 * Difference between both run modes is setting of SamplerBase::setRunLevel flag.
 */
void HerwigGenericRead(string reponame, string runname);

/**
 * Herwig read mode
 * 
 */
void HerwigRead(string reponame, string runname);

/**
 * Herwig build mode
 * 
 */
void HerwigBuild(string reponame, string runname);

/**
 * Herwig generic run mode
 * 
 * Function is used by integrate and run mode of Herwig++.
 * Difference between both run modes is either bool integrationJob is true or false.
 * Furthermore the SamplerBase::setRunLevel flag is set accordingly.
 */
void HerwigGenericRun(string runname, string setupfile,
	       int seed, string tag, long N, 
	       bool tics, bool resume, int jobs,
	       bool integrationJob,
	       string integrationList);

/**
 * Herwig run mode
 * 
 */
void HerwigRun(string runname, string setupfile,
	       int seed, string tag, long N, 
	       bool tics, bool resume, int jobs,
	       string integrationList);

/**
 * Herwig integration mode
 * 
 */
void HerwigIntegrate(string runname, string setupfile,
	       int seed, string tag, long N, 
	       bool tics, bool resume, int jobs,
	       string integrationList);

/**
 * Herwig++ main function
 * 
 * Event generation - it all starts from here. :-)
 * Function is only used in standalone mode, otherwise one can directly make use of the above run mode functions.
 */
int main(int argc, char * argv[]);

/**
 * Struct with enum to define runMode of Herwig++
 * 
 * Enum defines which runMode of Herwig++ should be used.
 */
struct HerwigRunMode {
  enum runMode { ERROR = -1, INIT = 1 , READ = 2, BUILD = 3, INTEGRATE = 4, RUN = 5 };
};

/**
 * Helper class to simplify read in and handling of command line parameters
 * 
 * Class is not needed if one links directly to the different Herwig++ run mode functions.
 * Class is defined as a singleton, so that setSearchPaths function can later be used in HerwigGenericRead function.
 */
class HelperReadInCommandLineParameters {
  private:
    /// Make default constructor private since argc and argv must be provided
    HelperReadInCommandLineParameters() {}
    
    /// Constructor is used but private since class should be a singleton.
    HelperReadInCommandLineParameters(int argc, char * argv[]);
    
    /// Forbid further instance creation by forbidding copy constructor
    HelperReadInCommandLineParameters( const HelperReadInCommandLineParameters& );
    
    /// Forbid further instance creation by copying instance
    HelperReadInCommandLineParameters& operator = (const HelperReadInCommandLineParameters &);
    
  public:
    
    ~HelperReadInCommandLineParameters() {cmdline_parser_free( &m_args_info );}
    
    /// Generates / returns singleton of this class
    static HelperReadInCommandLineParameters* instance(int argc, char * argv[])
    {
      if (!_instance)
	_instance = new HelperReadInCommandLineParameters(argc, argv);
      return _instance;
    }
    
    /*
     * Returns singleton after singleton was already created
     * 
     * Function is needed for use in HerwigRead to get a pre-initialized singleton instance
     */
    static HelperReadInCommandLineParameters* instance()
    {
      if (_instance && _instance->getReadInWasSuccessful()) {
	return _instance;
      } else {
	std::cerr << "Read in of command line parameters was not finished yet.";
	return NULL;
      }
    }
  
    /// Searches path of Herwig++ repo
    void setSearchPaths();
    
  private:
    /// static class object to make class a singleton, is initialized to null outside of class (see below)
    static HelperReadInCommandLineParameters* _instance;
    
    ///Gengetopt member variable which contains command line options
    gengetopt_args_info m_args_info;
       
    /// Interpret command status to define runMode
    int m_runMode;
    
    /// Check if read in of command line parameters was successful
    bool m_readInWasSuccessful;
    
  private:
    // Member variables
    /// Resume job
    bool m_resume;
    /// Activate tics
    bool m_tics;
    /// Name of in / run file 
    std::string m_runname;
    /// Name of repository
    std::string m_reponame;
    /// run name tag
    std::string m_tag;
    /// run modification filesystem
    std::string m_setupfile;
    /// List with integration jobs
    std::string m_integrationList;
    /// Number of events
    long m_N;
    /// Random number generator seed
    int m_seed;
    /// Number of parallel jobs if any
    int m_jobs;
    /// Number of integration jobs per cpu jobs
    unsigned int m_jobsize;
    /// Maximum number of integration jobs
    unsigned int m_maxjobs;
    
  public:		
    //Getter methods
    int getRunMode() {return m_runMode;}
    bool getReadInWasSuccessful() {return m_readInWasSuccessful;}
    bool getResume() {return m_resume;}
    bool getTics() {return m_tics;}
    std::string getRunName() {return m_runname;}
    std::string getRepoName() {return m_reponame;}
    std::string getTag() {return m_tag;}
    std::string getSetupFile() {return m_setupfile;}
    std::string getIntegrationList() {return m_integrationList;}
    long getN() {return m_N;}
    int getSeed() {return m_seed;}
    int getJobs() {return m_jobs;}
    unsigned int getJobSize() {return m_jobsize;}
    unsigned int getMaxJobs() {return m_maxjobs;} 
};
/// Initialize static singleton class object to Null
HelperReadInCommandLineParameters* HelperReadInCommandLineParameters::_instance = NULL;


#endif /* SRC_HERWIG_H */