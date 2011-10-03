// -*- C++ -*-
//
// SMFFPVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardModelFFPVertex class.
//

#include "SMFFPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

SMFFPVertex::SMFFPVertex()  : _charge(17,0.0), _couplast(0.), _q2last(0.*GeV2) {
  orderInGem(1);
  orderInGs(0);
}

void SMFFPVertex::doinit() {
  // PDG codes for the particles
  // the quarks
  for(int ix=1;ix<7;++ix) {
    addToList(-ix, ix, 22);
  }
  // the leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix, ix, 22);
  }
  for(int ix=1;ix<4;++ix) {
    _charge[2*ix-1]  = generator()->standardModel()->ed();
    _charge[2*ix ]   = generator()->standardModel()->eu();
    _charge[2*ix+9 ] = generator()->standardModel()->ee();
    _charge[2*ix+10] = generator()->standardModel()->enu();
  }
  FFVVertex::doinit();
}

void SMFFPVertex::persistentOutput(PersistentOStream & os) const {
  os << _charge;
}

void SMFFPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _charge;
}

ClassDescription<SMFFPVertex> 
SMFFPVertex::initSMFFPVertex;
// Definition of the static class description member.

void SMFFPVertex::Init() {
 static ClassDocumentation<SMFFPVertex> documentation
    ("The SMFFPVertex class is the implementation of"
     "the coupling of the photon to the Standard Model fermions");
}

// coupling for FFP vertex
void SMFFPVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr,tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = -electroMagneticCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  int iferm=abs(aa->id());
  if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
    left(_charge[iferm]);
    right(_charge[iferm]);
  }
  else throw HelicityConsistencyError() << "SMGFFPVertex::setCoupling "
					<< "Unknown particle in photon vertex" 
					<< Exception::runerror;
}
