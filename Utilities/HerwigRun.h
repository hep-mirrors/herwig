// -*- C++ -*-
#ifndef _HERWIG_RUN_H_
#define _HERWIG_RUN_H_

#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/StandardMatchers.h>
#include <ThePEG/PDT/DummyDecayer.h>
#include <ThePEG/Utilities/Debug.h>
#include "Herwig++/Utilities/HwDebug.h"
#include <ThePEG/Utilities/Timer.h>
#include <ThePEG/Utilities/DynamicLoader.h>
#include <ThePEG/Utilities/Exception.h>
#include <ThePEG/EventRecord/Event.h>
#include <ThePEG/Repository/Repository.h>
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
  enum RunStatus { UNKNOWN, INIT, READ, RUN };
  
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
   *  Destructor
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
   *  Names of files
   */
  //@{
  /**
   *  Name of the repository file
   */
  std::string repositoryFile() const;

  /**
   *  Name of the input file used to generate the repository
   */
  std::string repositoryInput() const;

  /**
   *   Name of the event generator to run
   */
  std::string runName() const;  
  //@}

  /**
   *  Access to various flags and parameters
   */
  //@{
  /**
   *  Number of events to generator
   */
  int getN() const;

  /**
   *  Number of events which have been generated
   */
  int getNGen() const;

  /**
   *  Status of the event generator
   */
  RunStatus status() const;

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

  /**
   *  Access to various particles and event properties
   */
  //@{
  /**
   *  Get the Step s from the event
   * @param e The event
   */
  ThePEG::StepVector getSteps(ThePEG::EventPtr e = ThePEG::EventPtr());

  /**
   *  Get the final-state particles from a Step
   * @param step The step, all steps if -1
   * @param e The event
   */
  ThePEG::tPVector getFinalState(int step = -1, ThePEG::EventPtr e = ThePEG::EventPtr());

  /**
   *  Get all the particle from a Step
   * @param step The step, all steps if -1
   * @param e The event
   */
  ThePEG::ParticleSet getAllParticles(int step = -1, 
				      ThePEG::EventPtr e = ThePEG::EventPtr());

  /**
   *  Get the intermediates from a Step
   * @param step The step, all steps if -1
   * @param e The event
   */
  ThePEG::ParticleSet getIntermediates(int step = -1, 
				       ThePEG::EventPtr e = ThePEG::EventPtr());

  /**
   *  Get the outgoing particles from a Step.
   * @param step The step, all steps if -1
   * @param e The event
   */
  ThePEG::ParticleSet getOutgoing(int step = -1, 
				  ThePEG::EventPtr e = ThePEG::EventPtr());
  //@}

private:

  /**
   *  default constructor is private to ensure arguements are passed
   */
  HerwigRun();

private:

  /**
   *  Number of events to generate
   */
  int N;

  /**
   *  Number of events which have been generated
   */
  int ngen;

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
   *  The timer for the run
   */
  ThePEG::MainTimer timer;

  /**
   *  Whether or not the event generator is initialised
   */
  bool isInitialized;

  /**
   *  Error status
   */
  bool errorFlag;

  /**
   *  The last event which was generated
   */
  ThePEG::EventPtr lastEvent;
};

}

#endif
