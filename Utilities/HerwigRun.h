#ifndef _HERWIG_RUN_H_
#define _HERWIG_RUN_H_

#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/StandardMatchers.h>
#include <ThePEG/PDT/PYDECYDummy.h>
#include <ThePEG/Utilities/Debug.h>
#include "Herwig++/Utilities/HwDebug.h"
#include <ThePEG/Utilities/Timer.h>
#include <ThePEG/Utilities/DynamicLoader.h>
#include <ThePEG/Utilities/Exception.h>
#include "Herwig++/Utilities/SmplHist.h"
#include <ThePEG/EventRecord/Event.h>
#include <ThePEG/Repository/Repository.h>
#include <iostream>

namespace Herwig {

class HerwigRun {
 public:
  enum RunStatus { UNKNOWN, INIT, READ, RUN };
 private:
  HerwigRun();
  int N;
  int ngen;
  int seed;
  std::string run;
  std::string repo;
  std::string repoin;
  bool egCreated;
  RunStatus Status;
  ThePEG::EGPtr eg;
  ThePEG::MainTimer timer;
  bool isInitialized;
  bool errorFlag;
  ThePEG::EventPtr lastEvent;

 public:
  HerwigRun(int argc, char **argv);
  ~HerwigRun();

  ThePEG::EGPtr eventGenerator();
  ThePEG::EventPtr generateEvent();
  std::string repositoryFile() const;
  std::string repositoryInput() const;
  std::string runName() const;  
  int getN() const;
  int getNGen() const;
  RunStatus status() const;
  bool isRunMode() const;
  bool isInitMode() const;
  bool isReadMode() const;
  static void printHelp(std::ostream &);
  ThePEG::StepVector getSteps(ThePEG::EventPtr e = ThePEG::EventPtr());
  ThePEG::tPVector getFinalState(int step = -1, ThePEG::EventPtr e = ThePEG::EventPtr());
  bool preparedToRun();
  ThePEG::ParticleSet getAllParticles(int step = -1, ThePEG::EventPtr e = ThePEG::EventPtr());
  ThePEG::ParticleSet getIntermediates(int step = -1, ThePEG::EventPtr e = ThePEG::EventPtr());
  ThePEG::ParticleSet getOutgoing(int step = -1, ThePEG::EventPtr e = ThePEG::EventPtr());
};

}

#endif
