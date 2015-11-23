// -*- C++ -*-
//
// EventShapesMasterAnalysis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EventShapesMasterAnalysis class.
//

#include "EventShapesMasterAnalysis.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void EventShapesMasterAnalysis::analyze(tEventPtr event, long ieve,
					int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
}

LorentzRotation EventShapesMasterAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void EventShapesMasterAnalysis::analyze(const tPVector & particles) {
  _shapes->reset(particles);
}

void EventShapesMasterAnalysis::analyze(tPPtr) {}

void EventShapesMasterAnalysis::persistentOutput(PersistentOStream & os) const {
  os << _shapes;
}

void EventShapesMasterAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> _shapes;
}

ClassDescription<EventShapesMasterAnalysis> EventShapesMasterAnalysis::initEventShapesMasterAnalysis;
// Definition of the static class description member.

void EventShapesMasterAnalysis::Init() {

  static ClassDocumentation<EventShapesMasterAnalysis> documentation
    ("The EventShapesMasterAnalysis class is the master class for event"
     " shapes analyses");

  static Reference<EventShapesMasterAnalysis,EventShapes> interfaceEventShapes
    ("EventShapes",
     "Pointer to the object which calculates the event shapes",
     &EventShapesMasterAnalysis::_shapes, false, false, true, false, false);

}

