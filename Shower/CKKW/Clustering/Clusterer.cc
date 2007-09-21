// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Clusterer class.
//

#include "Clusterer.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Clusterer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

Clusterer::~Clusterer() {}

void Clusterer::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _jetMeasure;
}

void Clusterer::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _jetMeasure;
}

AbstractClassDescription<Clusterer> Clusterer::initClusterer;
// Definition of the static class description member.

void Clusterer::Init() {

  static ClassDocumentation<Clusterer> documentation
    ("Clusterer defines the common interface for all clusterer objects "
     "used in reconstructing a parton shower history.");

}

