// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LEPEventShapes class.
//

#include "LEPEventShapes.h"
#include "Herwig++/Interfaces/EvtShapes.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LEPEventShapes.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LEPEventShapes::~LEPEventShapes() {}

void LEPEventShapes::analyze(tEventPtr event, long ieve, int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  // get the final-state particles
  tPVector hadrons=event->getFinalState();
  // event shapes
  Herwig::EvtShapes es(hadrons); 
  *_omthr += 1.-es.thrust();
  *_maj += es.thrustMajor();
  *_min += es.thrustMinor();
  *_obl += es.oblateness(); 
  *_c += es.CParameter(); 
  *_d += es.DParameter(); 
  *_sph += es.sphericity();
  *_apl += es.aplanarity();
  *_pla += es.planarity(); 
  *_mhi += es.Mhigh2();
  *_mlo += es.Mlow2(); 
  *_mdiff += es.Mdiff2(); 
  *_bmax += es.Bmax(); 
  *_bmin += es.Bmin(); 
  *_bsum += es.Bsum(); 
  *_bdiff += es.Bdiff(); 
}

LorentzRotation LEPEventShapes::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void LEPEventShapes::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void LEPEventShapes::analyze(tPPtr) {}

void LEPEventShapes::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void LEPEventShapes::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<LEPEventShapes> LEPEventShapes::initLEPEventShapes;
// Definition of the static class description member.

void LEPEventShapes::Init() {

  static ClassDocumentation<LEPEventShapes> documentation
    ("There is no documentation for the LEPEventShapes class");

}

