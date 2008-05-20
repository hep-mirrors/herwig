// -*- C++ -*-
//
// DecayConstructor.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayConstructor class.
//

#include "DecayConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr DecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr DecayConstructor::fullclone() const {
  return new_ptr(*this);
}
  
void DecayConstructor::persistentOutput(PersistentOStream & os) const {
   os << _theNBodyDecayConstructors;
}

void DecayConstructor::persistentInput(PersistentIStream & is, int) {
   is >> _theNBodyDecayConstructors;
}

ClassDescription<DecayConstructor> DecayConstructor::initDecayConstructor;
// Definition of the static class description member.

void DecayConstructor::Init() {

  static ClassDocumentation<DecayConstructor> documentation
    ("There is no documentation for the TwoBodyDecayConstructor class");

  static RefVector<DecayConstructor,Herwig::NBodyDecayConstructorBase> 
    interfaceNBodyDecayConstructors
    ("NBodyDecayConstructors",
     "Vector of references to N BodyDecayConstructors",
     &DecayConstructor::_theNBodyDecayConstructors, -1, false, false, true,
     false, false);

}

void DecayConstructor::createDecayers(const vector<PDPtr> & part) {
  for(unsigned int ix=0;ix < _theNBodyDecayConstructors.size();++ix) {
    _theNBodyDecayConstructors[ix]->init();
    _theNBodyDecayConstructors[ix]->DecayList(part);
  }
}

