// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LEPEventShapes class.
//

#include "LEPEventShapes.h"
#include "EvtShapes.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/Reference.h"
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
}

LorentzRotation LEPEventShapes::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void LEPEventShapes::analyze(const tPVector & particles) {
  *_omthr += 1.-_shapes->thrust();
  *_maj += _shapes->thrustMajor();
  *_min += _shapes->thrustMinor();
  *_obl += _shapes->oblateness(); 
  *_c += _shapes->CParameter(); 
  *_d += _shapes->DParameter(); 
  *_sph += _shapes->sphericity();
  *_apl += _shapes->aplanarity();
  *_pla += _shapes->planarity(); 
  *_mhi += _shapes->Mhigh2();
  *_mlo += _shapes->Mlow2(); 
  *_mdiff += _shapes->Mdiff2(); 
  *_bmax += _shapes->Bmax(); 
  *_bmin += _shapes->Bmin(); 
  *_bsum += _shapes->Bsum(); 
  *_bdiff += _shapes->Bdiff(); 
}

void LEPEventShapes::analyze(tPPtr) {}

void LEPEventShapes::persistentOutput(PersistentOStream & os) const {
  os << _shapes;
}

void LEPEventShapes::persistentInput(PersistentIStream & is, int) {
  is >> _shapes;
}

ClassDescription<LEPEventShapes> LEPEventShapes::initLEPEventShapes;
// Definition of the static class description member.

void LEPEventShapes::Init() {

  static ClassDocumentation<LEPEventShapes> documentation
    ("The LEPEventShapes class compares event shapes at the Z mass"
     "with experimental results");

  static Reference<LEPEventShapes,EventShapes> interfaceEventShapes
    ("EventShapes",
     "Pointer to the object which calculates the event shapes",
     &LEPEventShapes::_shapes, false, false, true, false, false);

}

