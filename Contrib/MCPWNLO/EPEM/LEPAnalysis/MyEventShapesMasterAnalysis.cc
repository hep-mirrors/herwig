// -*- C++ -*-
//
// MyEventShapesMasterAnalysis.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MyEventShapesMasterAnalysis class.
//

#include "MyEventShapesMasterAnalysis.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void MyEventShapesMasterAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
}

LorentzRotation MyEventShapesMasterAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void MyEventShapesMasterAnalysis::analyze(const tPVector & particles) {
  _shapes->reset(particles);
}

void MyEventShapesMasterAnalysis::analyze(tPPtr) {}

void MyEventShapesMasterAnalysis::persistentOutput(PersistentOStream & os) const {
  os << _shapes;
}

void MyEventShapesMasterAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> _shapes;
}

ClassDescription<MyEventShapesMasterAnalysis> MyEventShapesMasterAnalysis::initMyEventShapesMasterAnalysis;
// Definition of the static class description member.

void MyEventShapesMasterAnalysis::Init() {

  static ClassDocumentation<MyEventShapesMasterAnalysis> documentation
    ("The MyEventShapesMasterAnalysis class is the master class for event"
     " shapes analyses");

  static Reference<MyEventShapesMasterAnalysis,MyEventShapes> interfaceMyEventShapes
    ("MyEventShapes",
     "Pointer to the object which calculates the event shapes",
     &MyEventShapesMasterAnalysis::_shapes, false, false, true, false, false);

}

