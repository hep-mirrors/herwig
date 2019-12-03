// -*- C++ -*-
//
// RunDirectories.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_RunDirectories_H
#define HERWIG_RunDirectories_H
//
// This is the declaration of the RunDirectories class.
//

#include <string>
#include <list>

namespace Herwig {

/**
 * \author Simon Platzer
 *
 * \brief Handle directories for external library and grid storage.
 *
 */
class RunDirectories {

public:

  /**
   * Set a prefix for storing details of this run.
   */
  static void prefix(std::string p);

  /**
   * Return the prefix for storing details of this run.
   */
  static const std::string& prefix();

  /**
   * Return the name (and possibly create) a storage for build data
   */
  static const std::string& buildStorage();

  /**
   * Return true, if no run directories have been pushed yet
   */
  static bool empty();

  /**
   * Push a run identifier onto the run directories stack.
   */
  static void pushRunId(std::string);

  /**
   * Return (and possibly create) the top of the run directory stack
   * to be used for storage.
   */
  static const std::string& runStorage();

  /**
   * Return the storage to be used for interface order/contract files.
   */
  static const std::string& interfaceStorage();

public:

  /**
   * Default constructor fills the directory list to test.
   */
  RunDirectories();

  /**
   * Return true, if there are run directories still to be considered.
   */
  operator bool() const { return !directoriesLeft.empty(); }

  /**
   * Return true, if there are no run directories still to be considered.
   */
  bool operator!() const { return directoriesLeft.empty(); }

  /**
   * Return the next run directory to be considered and pop it from
   * the stack.
   */
  std::string nextRunStorage();

private:

  /**
   * The prefix for storing details of this run.
   */
  static std::string& thePrefix();

  /**
   * The build storage.
   */
  static std::string& theBuildStorage();

  /**
   * The list of run storage directories to be considered.
   */
  static std::list<std::string>& theRunDirectories();

  /**
   * The current run directory stack under consideration
   */
  std::list<std::string> directoriesLeft;

};

}

#endif /* HERWIG_RunDirectories_H */
