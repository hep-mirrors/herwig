// -*- C++ -*-
//
// SMFFZVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFZVertex class.
//

#include "SMFFZVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void SMFFZVertex::persistentOutput(PersistentOStream & os) const {
  os << _gl << _gr;
}

void SMFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _gl >> _gr;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SMFFZVertex,FFVVertex>
describeHerwigSMFFZVertex("Herwig::SMFFZVertex", "Herwig.so");

void SMFFZVertex::Init() {
  static ClassDocumentation<SMFFZVertex> documentation
    ("The SMFFZVertex class is the implementation of"
     "the coupling of the Z boson to the Standard Model fermions");
}

void SMFFZVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr,tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = -electroMagneticCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  int iferm=abs(aa->id());
  if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
    left(_gl[iferm]);
    right(_gr[iferm]);
  }
  else
    throw HelicityConsistencyError() << "SMFFZVertex::setCoupling "
				     << "Unknown particle in Z vertex" 
				     << Exception::runerror;
}

SMFFZVertex::SMFFZVertex() : _gl(17,0.0), _gr(17,0.0),
			     _couplast(0.0), _q2last(ZERO) {
  orderInGem(1);
  orderInGs(0);
}

void SMFFZVertex::doinit() {
  // PDG codes for the particles
  // the quarks
  for(int ix=1;ix<7;++ix) {
    addToList(-ix, ix, 23);
  }
  // the leptons
  for(int ix=11;ix<17;++ix) {
    addToList(-ix, ix, 23);
  }
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sin2ThetaW();
  double fact = 0.25/sqrt(sw2*(1.-sw2));
  for(int ix=1;ix<4;++ix) {
    _gl[2*ix-1]  = fact*(sm->vd()  + sm->ad() );
    _gl[2*ix ]   = fact*(sm->vu()  + sm->au() );
    _gl[2*ix+9 ] = fact*(sm->ve()  + sm->ae() );
    _gl[2*ix+10] = fact*(sm->vnu() + sm->anu());
    _gr[2*ix-1]  = fact*(sm->vd()  - sm->ad() );
    _gr[2*ix ]   = fact*(sm->vu()  - sm->au() );
    _gr[2*ix+9 ] = fact*(sm->ve()  - sm->ae() );
    _gr[2*ix+10] = fact*(sm->vnu() - sm->anu());
  }
  FFVVertex::doinit();
}
