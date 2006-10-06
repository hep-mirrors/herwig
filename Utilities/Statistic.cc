// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Statistic class.
//

#include "Statistic.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Statistic.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

Statistic::~Statistic() {}

void Statistic::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void Statistic::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<Statistic> Statistic::initStatistic;
// Definition of the static class description member.

void Statistic::Init() {

  static ClassDocumentation<Statistic> documentation
    ("There is no documentation for the Statistic class");

}

