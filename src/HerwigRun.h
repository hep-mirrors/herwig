// -*- C++ -*-
//
// HerwigRun.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_RUN_H_
#define HERWIG_RUN_H_

#include <ostream>

namespace Herwig {

/** \ingroup Utilities
 *
 *  This class is used to setup all the default command line arguments for
 *  Herwig++. This class provides a simple interface to many of ThePEG 
 *  structures used to store the event. 
 */
class HerwigRun {
 
public:

  /**
   *  Enumeration for the status of the event generator
   */  
  enum RunStatus { UNKNOWN, INIT, READ, RUN, ERROR, SUCCESS };
  
public:

  /**
   *   Constructors and Destructors
   */
  //@{
  /**
   *  Constructor with parameters
   * @param argc Length of string
   * @param argv string
   */
  HerwigRun(int argc, char **argv);

  /**
   * Destructor
   */
  ~HerwigRun();
  //@}
  
  /**
   * Generate the requested number of events
   */
  void generateEvents();

  /**
   *  Access to various flags and parameters
   */
  //@{
  /**
   *  State of the object
   */
  bool good() const;
  //@}

  /**
   *  Print help information on how to use the class
   */
  static void printHelp(std::ostream &);

private:

  /**
   *  Number of events to generate
   */
  long N_;

  /**
   *  Random number seed
   */
  int seed_;

  /**
   *  Name of event generator to run
   */
  std::string runname_;

  /**
   *  Status of the event generator
   */
  RunStatus status_;
};

}

#endif
