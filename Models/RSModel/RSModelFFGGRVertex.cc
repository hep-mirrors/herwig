// -*- C++ -*-
//
// RSModelFFGGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelFFGGRVertex class.
//

#include "RSModelFFGGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

RSModelFFGGRVertex::RSModelFFGGRVertex() 
  : couplast_(0.), q2last_(ZERO), kappa_(ZERO) {
  orderInGem(1);
  orderInGs (1);
}

void RSModelFFGGRVertex::doinit() {
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,ix,21,39);
  }
  FFVTVertex::doinit();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS) throw Exception() 
	      << "Must have RSModel in RSModelFFGGRVertex::doinit()"
	      << Exception::runerror;
  kappa_ = 2./hwRS->lambda_pi();
}

void RSModelFFGGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV);
}

void RSModelFFGGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV);
}

ClassDescription<RSModelFFGGRVertex> RSModelFFGGRVertex::initRSModelFFGGRVertex;
// Definition of the static class description member.

void RSModelFFGGRVertex::Init() {
  static ClassDocumentation<RSModelFFGGRVertex> documentation
    ("The RSModelFFGGRVertexxs class is the implementation"
     " of the two fermion vector coupling for the RS model.");
  
}

// FFGGR coupling
#ifndef NDEBUG
void RSModelFFGGRVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr,
				     tcPDPtr cc, tcPDPtr) {
#else
void RSModelFFGGRVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,
				      tcPDPtr, tcPDPtr) {
#endif
  // work out the particles
  assert(cc->id()==ParticleID::g && abs(aa->id()) <=6 );
  // overall factor
  if(q2last_ != q2 || couplast_ == 0. ) {
    couplast_ = strongCoupling(q2);
    q2last_ = q2;
  }
  left (1.);
  right(1.);
  // set the coupling
  norm( UnitRemoval::E * kappa_ * couplast_);
}
