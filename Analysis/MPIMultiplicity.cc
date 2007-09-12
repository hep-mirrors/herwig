// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MPIMultiplicity class.
//

#include "MPIMultiplicity.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MPIMultiplicity.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

MPIMultiplicity::~MPIMultiplicity() {}

void MPIMultiplicity::analyze(tEventPtr event, long , int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  if(!theShowerHandler->IsMPIOn()) return;
  theRealMult += generator()->eventHandler()->currentCollision()->subProcesses().size()-1;
  theRequestedMult += theMPIHandler->multiplicity();
}

LorentzRotation MPIMultiplicity::transform(tEventPtr ) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void MPIMultiplicity::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void MPIMultiplicity::analyze(tPPtr) {}

void MPIMultiplicity::persistentOutput(PersistentOStream & os) const {
  os << theShowerHandler << theMPIHandler;
}

void MPIMultiplicity::persistentInput(PersistentIStream & is, int) {
  is >> theShowerHandler >> theMPIHandler;
}

ClassDescription<MPIMultiplicity> MPIMultiplicity::initMPIMultiplicity;
// Definition of the static class description member.

void MPIMultiplicity::Init() {

  static ClassDocumentation<MPIMultiplicity> documentation
    ("There is no documentation for the MPIMultiplicity class");

  
  static Reference<MPIMultiplicity,ShowerHandler> interfaceShowerHandler
    ("ShowerHandler",
     "A reference to the ShowerHandler",
     &MPIMultiplicity::theShowerHandler, true, false, true, false, false);
  
  static Reference<MPIMultiplicity,MPIHandler> interfaceMPIHandler
    ("MPIHandler",
     "A reference to the MPIHandler",
     &MPIMultiplicity::theMPIHandler, true, false, true, false, false);

}

