// -*- C++ -*-
//
// ADDModelGGGGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelGGGGRVertex class.
//

#include "ADDModelGGGGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

ADDModelGGGGRVertex::ADDModelGGGGRVertex() 
  : kappa_(ZERO), r_(ZERO), couplast_(0.), q2last_(ZERO) {
  orderInGem(1);
  orderInGs (1);
  colourStructure(ColourStructure::SU3F);
}

void ADDModelGGGGRVertex::doinit() {
  addToList(21, 21, 21, 39);
  VVVTVertex::doinit();
  // set the graviton coupling 
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD) 
    throw Exception() 
      << "Must have ADDModel in ADDModelGGGGRVertex::doinit()"
      << Exception::runerror;
  kappa_=2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

void ADDModelGGGGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << ounit(r_,GeV);
}

void ADDModelGGGGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> iunit(r_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModelGGGGRVertex,VVVTVertex>
describeHerwigADDModelGGGGRVertex("Herwig::ADDModelGGGGRVertex", "HwADDModel.so");

void ADDModelGGGGRVertex::Init() {
 static ClassDocumentation<ADDModelGGGGRVertex> documentation
    ("The ADDModelGGGGRVertex class is the four point coupling"
     " of three vector bosons and a graviton in the Randell-Sundrum model.");
  
}

#ifndef NDEBUG
void ADDModelGGGGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
				     tcPDPtr c, tcPDPtr) {
#else
void ADDModelGGGGRVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,
				     tcPDPtr, tcPDPtr) {
#endif
  assert(a->id() == ParticleID::g && b->id() ==  ParticleID::g &&
	 c->id() == ParticleID::g);
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = strongCoupling(q2);
    q2last_=q2;
  }
  norm(Complex(couplast_*kappa_*UnitRemoval::E));
}

Complex ADDModelGGGGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
					Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}
