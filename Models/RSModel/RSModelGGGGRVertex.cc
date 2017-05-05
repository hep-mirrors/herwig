// -*- C++ -*-
//
// RSModelGGGGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelGGGGRVertex class.
//

#include "RSModelGGGGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

RSModelGGGGRVertex::RSModelGGGGRVertex() 
  : kappa_(ZERO), _couplast(0.), _q2last(ZERO) {
  orderInGem(1);
  orderInGs (1);
}

void RSModelGGGGRVertex::doinit() {
  addToList(21, 21, 21, 39);
  VVVTVertex::doinit();
  // set the graviton coupling 
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS) 
    throw Exception() 
      << "Must have RSModel in RSModelGGGGRVertex::doinit()"
      << Exception::runerror;
  kappa_=2./hwRS->lambda_pi();
}

void RSModelGGGGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV);
}
void RSModelGGGGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV);
}

ClassDescription<RSModelGGGGRVertex> RSModelGGGGRVertex::initRSModelGGGGRVertex;
// Definition of the static class description member.

void RSModelGGGGRVertex::Init() {
 static ClassDocumentation<RSModelGGGGRVertex> documentation
    ("The RSModelGGGGRVertex class is the four point coupling"
     " of three vector bosons and a graviton in the Randell-Sundrum model.");
  
}


// couplings for the GGGGR vertex
#ifndef NDEBUG
void RSModelGGGGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
				     tcPDPtr c, tcPDPtr) {
#else
void RSModelGGGGRVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,
				     tcPDPtr, tcPDPtr) {
#endif
  assert(a->id() == ParticleID::g && b->id() ==  ParticleID::g &&
	 c->id() == ParticleID::g);
  if(q2!=_q2last||_couplast==0.) {
    _couplast = strongCoupling(q2);
    _q2last=q2;
  }
  norm(Complex(_couplast*kappa_*UnitRemoval::E));
}
