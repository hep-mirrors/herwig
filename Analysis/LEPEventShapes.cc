// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LEPEventShapes class.
//

#include "LEPEventShapes.h"
#include "EventShapes.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LEPEventShapes.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LEPEventShapes::analyze(tEventPtr event, long ieve, int loop, int state) {
  eventweight_ = event->weight();
  AnalysisHandler::analyze(event, ieve, loop, state);
  if ( loop > 0 || state != 0 || !event ) return;
  // get the final-state particles
  tPVector hadrons=event->getFinalState();
  // event shapes
}

LorentzRotation LEPEventShapes::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void LEPEventShapes::analyze(const tPVector & particles) {
  _shapes->reset(particles);
  _omthr ->addWeighted( 1.-_shapes->thrust() ,eventweight_);
  _maj ->addWeighted( _shapes->thrustMajor() ,eventweight_);
  _min ->addWeighted( _shapes->thrustMinor() ,eventweight_);
  _obl ->addWeighted( _shapes->oblateness() ,eventweight_); 
  _c ->addWeighted( _shapes->CParameter() ,eventweight_); 
  _d ->addWeighted( _shapes->DParameter() ,eventweight_); 
  _sph ->addWeighted( _shapes->sphericity() ,eventweight_);
  _apl ->addWeighted( _shapes->aplanarity() ,eventweight_);
  _pla ->addWeighted( _shapes->planarity() ,eventweight_); 
  _mhi ->addWeighted( _shapes->Mhigh2() ,eventweight_);
  _mlo ->addWeighted( _shapes->Mlow2() ,eventweight_); 
  _mdiff ->addWeighted( _shapes->Mdiff2() ,eventweight_); 
  _bmax ->addWeighted( _shapes->Bmax() ,eventweight_); 
  _bmin ->addWeighted( _shapes->Bmin() ,eventweight_); 
  _bsum ->addWeighted( _shapes->Bsum() ,eventweight_); 
  _bdiff ->addWeighted( _shapes->Bdiff() ,eventweight_); 
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

