// -*- C++ -*-
//
// SMFFZVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFZVertex class.
//

#include "SMFFZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void SMFFZVertex::persistentOutput(PersistentOStream & os) const {
  os << _gl << _gr <<  _theSM;
}

void SMFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _gl >> _gr >> _theSM;
}

ClassDescription<SMFFZVertex> 
SMFFZVertex::initSMFFZVertex;
// Definition of the static class description member.

void SMFFZVertex::Init() {
  static ClassDocumentation<SMFFZVertex> documentation
    ("The SMFFZVertex class is the implementation of"
     "the coupling of the Z boson to the Standard Model fermions");
}

void SMFFZVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr,tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last) {
    double alpha = _theSM->alphaEM(q2);
    _couplast = -sqrt(4.0*Constants::pi*alpha);
    _q2last=q2;
  }
  setNorm(_couplast);
  // the left and right couplings
  int iferm=abs(a->id());
  if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
    setLeft(_gl[iferm]);
    setRight(_gr[iferm]);
  }
  else
    throw HelicityConsistencyError() << "SMFFZVertex::setCoupling "
				     << "Unknown particle in Z vertex" 
				     << Exception::runerror;
}

SMFFZVertex::SMFFZVertex() : _gl(17,0.0), _gr(17,0.0),
			     _couplast(0.0), _q2last(0.*GeV2) {
  // PDG codes for the particles
  vector<int> first,second,third;
  // the quarks
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(23);
  }
  // the leptons
  for(unsigned int ix=11;ix<17;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(23);
  }
  setList(first,second,third);
}

void SMFFZVertex::doinit() throw(InitException) {
  _theSM = generator()->standardModel();
  double sw2=_theSM->sin2ThetaW();
  double fact = 0.25/sqrt(sw2*(1.-sw2));
  for(int ix=1;ix<4;++ix) {
    _gl[2*ix-1]  = fact*(_theSM->vd()  + _theSM->ad() );
    _gl[2*ix ]   = fact*(_theSM->vu()  + _theSM->au() );
    _gl[2*ix+9 ] = fact*(_theSM->ve()  + _theSM->ae() );
    _gl[2*ix+10] = fact*(_theSM->vnu() + _theSM->anu());
    _gr[2*ix-1]  = fact*(_theSM->vd()  - _theSM->ad() );
    _gr[2*ix ]   = fact*(_theSM->vu()  - _theSM->au() );
    _gr[2*ix+9 ] = fact*(_theSM->ve()  - _theSM->ae() );
    _gr[2*ix+10] = fact*(_theSM->vnu() - _theSM->anu());
  }
  orderInGem(1);
  orderInGs(0);
  FFVVertex::doinit();
}
