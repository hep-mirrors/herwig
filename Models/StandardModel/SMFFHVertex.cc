// -*- C++ -*-
//
// SMFFHVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFHVertex class.
//

#include "SMFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;

void SMFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << _sw;
}

void SMFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV) >> _sw;
}

ClassDescription<SMFFHVertex> 
SMFFHVertex::initSMFFHVertex;
// Definition of the static class description member.

void SMFFHVertex::Init() {

  static ClassDocumentation<SMFFHVertex> documentation
    ("The SMFFHVertex class is the implementation"
     " of the helicity amplitude calculation of the Standard Model Higgs"
     " fermion-antiferiom vertex.");
  
}

void SMFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr, tcPDPtr, int) {
  int iferm=abs(a->id());
  // left and right couplings set to one
  setLeft(1.); setRight(1.);
  // first the overall normalisation
  if(q2!=_q2last) {
    double alpha = _theSM->alphaEM(q2);
    _couplast = -0.5*sqrt(4.0*Constants::pi*alpha)/_sw/_mw;
    _q2last=q2;
    _idlast=iferm;
    if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
      _masslast=_theSM->mass(q2,a);
    }
    else {
      throw HelicityConsistencyError() << "SMFFHVertex::setCoupling " 
				       << "Unknown particle in Higgs vertex" 
				       << Exception::warning;
      _masslast = 0*MeV;
    }
  }
  else if(iferm!=_idlast) {
    _idlast=iferm;
    if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
      _masslast=_theSM->mass(q2,a);
    }
    else {
      throw HelicityConsistencyError() << "SMFFHVertex::setCoupling " 
				       << "Unknown particle in Higgs vertex" 
				       << Exception::warning;
      _masslast = 0*MeV;
    }
  }
  setNorm(_couplast*_masslast);
}
