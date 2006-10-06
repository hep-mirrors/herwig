// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonExtractor class.
//

#include "PartonExtractor.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PartonExtractor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

Herwig::PartonExtractor::~PartonExtractor() {}

void Herwig::PartonExtractor::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void Herwig::PartonExtractor::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<Herwig::PartonExtractor> Herwig::PartonExtractor::initPartonExtractor;
// Definition of the static class description member.

void Herwig::PartonExtractor::Init() {

  static ClassDocumentation<PartonExtractor> documentation
    ("There is no documentation for the PartonExtractor class");

}

void Herwig::PartonExtractor::
colourConnect(tPPtr particle, tPPtr parton, const tPVector & remnants) const {

  // Sorry cannot handle coloured resolved particles.
  if ( particle->coloured() ) throw RemColException(*this);

  // First connect the loose colour line from the extacted parton.
  if ( parton->hasColour() ) parton->colourLine()->addColoured(remnants[0],true);

  // First connect the loose anti-colour line from the extacted parton.
  if ( parton->hasAntiColour() )
    parton->antiColourLine()->addColoured(remnants[0],false);
}
