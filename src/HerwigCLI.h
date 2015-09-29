// -*- C++ -*-
//
// Herwig.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2015 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef SRC_HERWIG_CLI_H
#define SRC_HERWIG_CLI_H

#include "HerwigUI.h"

namespace Herwig {

class HerwigCLI : public HerwigUI {
public:
  HerwigCLI(int argc, char * argv[]);

  ~HerwigCLI();

  RunMode::Mode runMode() const { return runMode_; }

  bool resume() const { return resume_; }
  bool tics() const { return tics_; }
  std::string tag() const { return tag_; }

  std::string inputfile() const { return inputfile_; }
  std::string repository() const { return repository_; }
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

  void quitWithError() const;

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