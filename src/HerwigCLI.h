// -*- C++ -*-
//
// HerwigCLI.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2016 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef SRC_HERWIG_CLI_H
#define SRC_HERWIG_CLI_H

#include "Herwig/API/HerwigUI.h"
#include <iostream>

namespace Herwig {

/**
 * HerwigCLI is the default implementation of the HerwigUI interface.
 * 
 * Using the gengetopt tool, we fill all required pieces from reading
 * command line flags from the main executable.
 */
class HerwigCLI : public HerwigUI {
public:

  /// Constructor from the arguments provided by main()
  HerwigCLI(int argc, char * argv[]);

  /// Destructor to leave a clean ThePEG::Repository behind
  ~HerwigCLI();

  /// Requested Herwig run mode
  RunMode::Mode runMode() const { return runMode_; }

  /// Try to resume execution from an earlier interrupted run.
  bool resume() const { return resume_; }

  /// Require verbose progress markers
  bool tics() const { return tics_; }

  /// A user-defined tag to append to the run name.
  std::string tag() const { return tag_; }

  /// Name of the file to be read
  std::string inputfile() const { return inputfile_; }

  /// Repository name to operate on
  std::string repository() const { return repository_; }

  /// Name of the setup file to be read, to modify the repository
  std::string setupfile() const { return setupfile_; }
 
  std::string integrationList() const { return integrationList_; }


  const std::vector<std::string> & 
  prependReadDirectories() const { return prependReadDirectories_; }

  const std::vector<std::string> & 
  appendReadDirectories() const { return appendReadDirectories_; }

  long N() const { return N_; }
  int seed() const { return seed_; }
  int jobs() const { return jobs_; }
  unsigned int jobSize() const { return jobsize_; }
  unsigned int maxJobs() const { return maxjobs_; }  

  void quitWithHelp() const;

  void quit() const;

   /// Return the standard out stream to be used
  virtual std::ostream& outStream() const { return std::cout; }

   /// Return the standard err stream to be used
  virtual std::ostream& errStream() const { return std::cerr; }

  /// Return the standard in stream to be used
  virtual std::istream& inStream() const { return std::cin; }

private:

  RunMode::Mode runMode_;

  bool resume_;
  bool tics_;
  std::string tag_;

  std::string inputfile_;
  std::string repository_;
  std::string setupfile_;

  std::string integrationList_;

  std::vector<std::string> prependReadDirectories_;
  std::vector<std::string> appendReadDirectories_;

  long N_;
  int seed_;
  int jobs_;
  unsigned int jobsize_;
  unsigned int maxjobs_;

};

}

#endif
