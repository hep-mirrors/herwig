#ifndef _HERWIG_RUN_H_
#define _HERWIG_RUN_H_

#include "Pythia7/Repository/EventGenerator.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/PDT/StandardMatchers.h"
#include "Pythia7/PDT/PYDECYDummy.h"
#include "Pythia7/Utilities/Debug.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Utilities/Timer.h"
#include "Pythia7/Utilities/DynamicLoader.h"
#include "Pythia7/Misc/Exception.h"
#include "Herwig++/Utilities/SmplHist.h"
#include "Pythia7/EventRecord/Event.h"
#include "Pythia7/Repository/Repository.h" 
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
  Pythia7::EGPtr eg;
  Pythia7::MainTimer timer;
  bool isInitialized;
  bool errorFlag;
  Pythia7::EventPtr lastEvent;

 public:
  HerwigRun(int argc, char **argv);
  ~HerwigRun();

  Pythia7::EGPtr eventGenerator();
  Pythia7::EventPtr generateEvent();
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
  Pythia7::tPVector getFinalState(int step = -1, 
				  Pythia7::EventPtr e = Pythia7::EventPtr());
  bool preparedToRun();
};

}

#endif
