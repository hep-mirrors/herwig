// -*- C++ -*-
//
// Herwig.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef SRC_HERWIG_H
#define SRC_HERWIG_H

namespace Herwig {

class HerwigUI;

/**
 * A very high-level API for interacting with Herwig's different run modes.
 *
 * It's not very convenient (yet), since you'll have to provide your own 
 * HerwigUI-derived object with some fairly obscure settings.
 *
 * Much more fine-grained control is available through ThePEG::Repository.
 */
namespace API {

void init(const HerwigUI &);

void read(const HerwigUI &);

void build(const HerwigUI &);

void integrate(const HerwigUI &);

void run(const HerwigUI &);

}


}























/**
 * Show usage information in terminal
 * 
 */
//void printUsageAndExit();

/**
 * Searches path of Herwig repo
 * 
 * You can define two string vectors with directories which Herwig will use to look in for files.
 * A vector with directories which will be prepended and a vector with directories which will be appended.
 * By default both vectors are optional.
 * 
 * @param[in] prependReadDirectories Directories from which Herwig read files, directories will be prepended to read path
 * @param[in] appendReadDirectories Directories from which Herwig read files, directories will be appended to read path
 * @param[in] useCWD If current working directory should be prepended to search path.
 */
//void setSearchPaths(const std::vector<std::string> & prependReadDirectories = std::vector<std::string>(), 
//        const std::vector<std::string> & appendReadDirectories = std::vector<std::string>(),
//                    bool useCWD = true);

/**
 * Herwig init run mode
 * 
 * @param[in] infile Name of infile
 * @param[in] reponame Name of repository
 */
//void HerwigInit(string infile, string reponame);

/**
 * Herwig generic read mode
 * 
 * Function is used by built and read mode of Herwig.
 * Difference between both run modes is setting of SamplerBase::setRunLevel flag.
 * You can define two string vectors with directories which Herwig will use to look in for files.
 * A vector with directories which will be prepended and a vector with directories which will be appended.
 * By default both vectors are optional.
 * 
 * @param[in] reponame Name of repository
 * @param[in] runname Name of the run 
 * @param[in] prependReadDirectories Directories from which Herwig read files, directories will be prepended to read path
 * @param[in] appendReadDirectories Directories from which Herwig read files, directories will be appended to read path
 */
//void HerwigGenericRead(string reponame, string runname,
//    const std::vector<std::string> & prependReadDirectories = std::vector<std::string>(), 
//    const std::vector<std::string> & appendReadDirectories = std::vector<std::string>());

/**
 * Herwig read mode
 * 
 * You can define two string vectors with directories which Herwig will use to look in for files.
 * A vector with directories which will be prepended and a vector with directories which will be appended.
 * By default both vectors are optional.
 * 
 * @param[in] reponame Name of repository
 * @param[in] runname Name of the run 
 * @param[in] prependReadDirectories Directories from which Herwig read files, directories will be prepended to read path
 * @param[in] appendReadDirectories Directories from which Herwig read files, directories will be appended to read path
 */
//void HerwigRead(string reponame, string runname,
//    const std::vector<std::string> & prependReadDirectories = std::vector<std::string>(), 
//    const std::vector<std::string> & appendReadDirectories = std::vector<std::string>());

/**
 * Herwig build mode
 *
 * You can define two string vectors with directories which Herwig will use to look in for files.
 * A vector with directories which will be prepended and a vector with directories which will be appended.
 * By default both vectors are optional. 
 * 
 * @param[in] reponame Name of repository
 * @param[in] runname Name of the run
 * @param[in] prependReadDirectories Directories from which Herwig read files, directories will be prepended to read path
 * @param[in] appendReadDirectories Directories from which Herwig read files, directories will be appended to read path
 */
//void HerwigBuild(string reponame, string runname,
//    const std::vector<std::string> & prependReadDirectories = std::vector<std::string>(), 
//    const std::vector<std::string> & appendReadDirectories = std::vector<std::string>());

/**
 * Herwig generic run mode
 * 
 * Function is used by integrate and run mode of Herwig.
 * Difference between both run modes is either bool integrationJob is true or false.
 * Furthermore the SamplerBase::setRunLevel flag is set accordingly.
 * 
 * @param[in] runname Name of the run
 * @param[in] setupfile Name of an optional setupfile
 * @param[in] seed Seed value for the random number generator
 * @param[in] tag Name of the run tag
 * @param[in] N Number of events to generate
 * @param[in] tics Defines if tics should be shown
 * @param[in] resume Defines if previous job should be resumed
 * @param[in] jobs Defines how many cpu jobs should be used
 * @param[in] integrationJob Defines if either integration or run runMode should be activated
 * @param[in] integrationList List containing the integration jobs 
 */
//void HerwigGenericRun(string runname, string setupfile,
 //        int seed, string tag, long N, 
//         bool tics, bool resume, int jobs,
 //        bool integrationJob,
 //        string integrationList);

/**
 * Herwig run mode
 * 
 * @param[in] runname Name of the run
 * @param[in] setupfile Name of an optional setupfile
 * @param[in] seed Seed value for the random number generator
 * @param[in] tag Name of the run tag
 * @param[in] N Number of events to generate
 * @param[in] tics Defines if tics should be shown
 * @param[in] resume Defines if previous job should be resumed
 * @param[in] jobs Defines how many cpu jobs should be used
 * @param[in] integrationList List containing the integration jobs 
 */
//void HerwigRun(string runname, string setupfile,
//         int seed, string tag, long N, 
//         bool tics, bool resume, int jobs,
//         string integrationList);

/**
 * Herwig integration mode
 * 
 * @param[in] runname Name of the run
 * @param[in] setupfile Name of an optional setupfile
 * @param[in] seed Seed value for the random number generator
 * @param[in] tag Name of the run tag
 * @param[in] N Number of events to generate
 * @param[in] tics Defines if tics should be shown
 * @param[in] resume Defines if previous job should be resumed
 * @param[in] jobs Defines how many cpu jobs should be used
 * @param[in] integrationJob Defines if either integration or run runMode should be activated
 * @param[in] integrationList List containing the integration jobs 
 */
//void HerwigIntegrate(string runname, string setupfile,
//         int seed, string tag, long N, 
//         bool tics, bool resume, int jobs,
//         string integrationList);


#endif
