// -*- C++ -*-
//
// RunDirectories.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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

using std::string;
using std::list;

/**
 * \ingroup Matchbox
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
  static void prefix(string p);

  /**
   * Return the prefix for storing details of this run.
   */
  static const string& prefix();

  /**
   * Return the name (and possibly create) a storage for build data
   */
  static const string& buildStorage();

  /**
   * Return true, if no run directories have been pushed yet
   */
  static bool empty();

  /**
   * Push a run identifier onto the run directories stack.
   */
  static void pushRunId(string);

  /**
   * Return (and possibly create) the top of the run directory stack
   * to be used for storage.
   */
  static const string& runStorage();

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
  string nextRunStorage();

private:

  /**
   * The prefix for storing details of this run.
   */
  static string& thePrefix();

  /**
   * The build storage.
   */
  static string& theBuildStorage();

  /**
   * The list of run storage directories to be considered.
   */
  static list<string>& theRunDirectories();

  /**
   * The current run directory stack under consideration
   */
  list<string> directoriesLeft;

};

}

#endif /* HERWIG_RunDirectories_H */
