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
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

void LEPEventShapes::analyze(tEventPtr event, long ieve, int loop, int state) {
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
  double wgt=CurrentGenerator::current().currentEvent()->weight();
  _shapes->reset(particles);
  _omthr ->addWeighted( 1.-_shapes->thrust() ,wgt);
  _maj ->addWeighted( _shapes->thrustMajor() ,wgt);
  _min ->addWeighted( _shapes->thrustMinor() ,wgt);
  _obl ->addWeighted( _shapes->oblateness() ,wgt); 
  _c ->addWeighted( _shapes->CParameter() ,wgt); 
  _d ->addWeighted( _shapes->DParameter() ,wgt); 
  _sph ->addWeighted( _shapes->sphericity() ,wgt);
  _apl ->addWeighted( _shapes->aplanarity() ,wgt);
  _pla ->addWeighted( _shapes->planarity() ,wgt); 
  _mhi ->addWeighted( _shapes->Mhigh2() ,wgt);
  _mlo ->addWeighted( _shapes->Mlow2() ,wgt); 
  _mdiff ->addWeighted( _shapes->Mdiff2() ,wgt); 
  _bmax ->addWeighted( _shapes->Bmax() ,wgt); 
  _bmin ->addWeighted( _shapes->Bmin() ,wgt); 
  _bsum ->addWeighted( _shapes->Bsum() ,wgt); 
  _bdiff ->addWeighted( _shapes->Bdiff() ,wgt); 
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

