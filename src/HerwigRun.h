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

#include <ThePEG/EventRecord/Step.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <iostream>

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
   *  Access to the event generator pointer
   */
  ThePEG::EGPtr eventGenerator();

  /**
   *  generate an event
   */
  ThePEG::EventPtr generateEvent();

  /**
   * Generate the requested number of events
   */
  void generateEvents();

  /**
   *  Access to various flags and parameters
   */
  //@{
  /**
   *  Number of events to generator
   */
  long getN() const;

  /**
   *  Number of events which have been generated
   */
  long getNGen() const;

  /**
   *  Status of the event generator
   */
  RunStatus status() const;

  /**
   *  State of the object
   */
  bool good() const;

  /**
   *  Is this the run mode
   */
  bool isRunMode() const;

  /**
   *  Is this the init mode
   */
  bool isInitMode() const;

  /**
   *  Is this the read mode
   */
  bool isReadMode() const;

  /**
   *  Is the event generator ready to run
   */
  bool preparedToRun();
  //@}

  /**
   *  Print help information on how to use the class
   */
  static void printHelp(std::ostream &);

private:

  /**
   *  Number of events to generate
   */
  long N;

  /**
   *  Number of events which have been generated
   */
  long ngen;

  /**
   *  Random number seed
   */
  int seed;

  /**
   *  Name of event generator to run
   */
  std::string run;

  /**
   *  Name of the repository file
   */
  std::string repo;

  /**
   *  Name of the input file to generate the repository
   */
  std::string repoin;

  /**
   *  Whether or not the event generator has been created
   */
  bool egCreated;

  /**
   *  Status of the event generator
   */
  RunStatus Status;

  /**
   *  Pointer to the event generator
   */
  ThePEG::EGPtr eg;

  /**
   *  Whether or not the event generator is initialised
   */
  bool isInitialized;

  /**
   *  The last event which was generated
   */
  ThePEG::EventPtr lastEvent;
};

}

#endif
