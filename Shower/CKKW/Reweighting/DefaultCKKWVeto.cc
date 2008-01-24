// -*- C++ -*-
//
// DefaultCKKWVeto.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DefaultCKKWVeto class.
//

#include "DefaultCKKWVeto.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DefaultCKKWVeto.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DefaultCKKWVeto::~DefaultCKKWVeto() {}

void DefaultCKKWVeto::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _resolution << _enabled << _optResolution;
}

void DefaultCKKWVeto::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _resolution >> _enabled >> _optResolution;
}

ClassDescription<DefaultCKKWVeto> DefaultCKKWVeto::initDefaultCKKWVeto;
// Definition of the static class description member.

void DefaultCKKWVeto::Init() {

  static ClassDocumentation<DefaultCKKWVeto> documentation
    ("DefaultCKKWVeto is a ShowerVeto wrapper around DefaultJetMeasure.");

}

