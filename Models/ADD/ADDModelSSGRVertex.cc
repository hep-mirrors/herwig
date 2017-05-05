// -*- C++ -*-
//
// ADDModelSSGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelSSGRVertex class.
//

#include "ADDModelSSGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

ADDModelSSGRVertex::ADDModelSSGRVertex() : kappa_(ZERO), r_(ZERO) {
  orderInGem(1);
  orderInGs (0);
}

void ADDModelSSGRVertex::doinit() {
  addToList(25,25,39);
  SSTVertex::doinit();
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD) 
    throw Exception() << "Must have ADDModel in ADDModelSSGRVertex::doinit()"
		      << Exception::runerror;
  kappa_=2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

void ADDModelSSGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << ounit(r_,GeV);
}

void ADDModelSSGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> iunit(r_,GeV);
}

ClassDescription<ADDModelSSGRVertex> ADDModelSSGRVertex::initADDModelSSGRVertex;
// Definition of the static class description member.

void ADDModelSSGRVertex::Init() {
  static ClassDocumentation<ADDModelSSGRVertex> documentation
    ("The ADDModelSSGRVertex class is the implementation of"
     " the ADDModel scalar-scalar-graviton vertex");
  
}

void ADDModelSSGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
  norm(Complex(kappa_ * UnitRemoval::E));
}

Complex ADDModelSSGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
				       Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}
