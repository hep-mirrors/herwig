// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DefaultJetMeasure class.
//

#include "DefaultJetMeasure.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DefaultJetMeasure.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DefaultJetMeasure::~DefaultJetMeasure() {}

void DefaultJetMeasure::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void DefaultJetMeasure::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

AbstractClassDescription<DefaultJetMeasure> DefaultJetMeasure::initDefaultJetMeasure;
// Definition of the static class description member.

void DefaultJetMeasure::Init() {

  static ClassDocumentation<DefaultJetMeasure> documentation
    ("DefaultJetMeasure is a jet resolution for the standard "
     "CKKW approaches together with the standard shower implementation. "
     "It provides methods to veto shower emissions, returns phase "
     "space boundaries as well as Jacobians for Sudakov exponent "
     "integration.");

}

