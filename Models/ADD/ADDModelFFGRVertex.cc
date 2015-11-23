// -*- C++ -*-
//
// ADDModelFFGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelFFGRVertex class.
//

#include "ADDModelFFGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void ADDModelFFGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << ounit(r_,GeV);
}

void ADDModelFFGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> iunit(r_,GeV);
}

ClassDescription<ADDModelFFGRVertex> ADDModelFFGRVertex::initADDModelFFGRVertex;
// Definition of the static class description member.

void ADDModelFFGRVertex::Init() {
  static ClassDocumentation<ADDModelFFGRVertex> documentation
    ("The ADDModelFFGRVertex class is the ADDModel calculation"
     " of the fermion-antifermion-graviton vertex");
  
}
  
void ADDModelFFGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
  norm(Complex(kappa_ * UnitRemoval::E));
}

ADDModelFFGRVertex::ADDModelFFGRVertex() : kappa_(ZERO), r_(ZERO) {
  orderInGem(1);
  orderInGs (0);
}

void ADDModelFFGRVertex::doinit() {
  // PDG codes for the particles
  // the quarks
  for (int ix=1;ix<7;++ix) addToList(-ix,ix,39);
  // the leptons
  for (int ix=11;ix<17;++ix) addToList(-ix,ix,39);
  FFTVertex::doinit();
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD)
    throw Exception() << "Must have ADDModel in ADDModelFFGRVertex::doinit()"
		      << Exception::runerror;
  kappa_=2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

Complex ADDModelFFGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
				       Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}
