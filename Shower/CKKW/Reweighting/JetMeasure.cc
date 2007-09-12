// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the JetMeasure class.
//

#include "JetMeasure.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "JetMeasure.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

JetMeasure::~JetMeasure() {}

void JetMeasure::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << ounit(_resolution,GeV2) << _min << _max;
}

void JetMeasure::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> iunit(_resolution,GeV2) >> _min >> _max;
}

AbstractClassDescription<JetMeasure> JetMeasure::initJetMeasure;
// Definition of the static class description member.

void JetMeasure::Init() {

  static ClassDocumentation<JetMeasure> documentation
    ("JetMeasure is the base class for evaluating Sudakovs, applying ME cuts"
     "and vetoing shower emissions in a CKKW ME+PS approach.");


  static Parameter<JetMeasure,Energy2> interfaceResolution
    ("Resolution",
     "The matching scale",
     &JetMeasure::_resolution, GeV2, 100.0*GeV2, 0*GeV2, 0*GeV2,
     false, false, Interface::nolimits);


}

