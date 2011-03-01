// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#include "EventExtractor.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EventExtractor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Analysis2;

EventExtractor::~EventExtractor() {}

void EventExtractor::persistentOutput(PersistentOStream & ) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void EventExtractor::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<EventExtractor> EventExtractor::initEventExtractor;
// Definition of the static class description member.

void EventExtractor::Init() {

  static ClassDocumentation<EventExtractor> documentation
    ("EventExtractor");


}


void EventExtractor::prepare() {

  _lastWeight = lastEvent()->weight();

  _lastMomenta.clear();

  tPVector parts = lastEvent()->getFinalState();

  for (tPVector::const_iterator p = parts.begin();
       p != parts.end(); ++p)
    _lastMomenta.push_back((**p).momentum());

  _ok = true;
  _doneEvent = false;

}
