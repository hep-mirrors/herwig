// -*- C++ -*-
//
// ADDModelFFGGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelFFGGRVertex class.
//

#include "ADDModelFFGGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


using namespace Herwig;
using namespace ThePEG;

ADDModelFFGGRVertex::ADDModelFFGGRVertex() 
  : couplast_(0.), q2last_(ZERO), kappa_(ZERO), r_(ZERO) {
  orderInGem(1);
  orderInGs (1);
  colourStructure(ColourStructure::SU3TFUND);
}

void ADDModelFFGGRVertex::doinit() {
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,ix,21,39);
  }
  FFVTVertex::doinit();
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD) throw Exception() 
	      << "Must have ADDModel in ADDModelFFGGRVertex::doinit()"
	      << Exception::runerror;
  kappa_ = 2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

void ADDModelFFGGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << ounit(r_,GeV);
}

void ADDModelFFGGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> iunit(r_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModelFFGGRVertex,FFVTVertex>
describeHerwigADDModelFFGGRVertex("Herwig::ADDModelFFGGRVertex", "HwADDModel.so");

void ADDModelFFGGRVertex::Init() {
  static ClassDocumentation<ADDModelFFGGRVertex> documentation
    ("The ADDModelFFGGRVertexxs class is the implementation"
     " of the two fermion vector coupling for the ADD model.");
  
}

#ifndef NDEBUG
void ADDModelFFGGRVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr,
				      tcPDPtr cc, tcPDPtr) {
#else
void ADDModelFFGGRVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,
				      tcPDPtr, tcPDPtr) {
#endif
  // work out the particles
  assert(cc->id()==ParticleID::g && abs(aa->id()) <= 6);
  // overall factor
  if(q2last_!=q2||couplast_==0.) {
    couplast_ = strongCoupling(q2);
    q2last_=q2;
  }
  left (1.);
  right(1.);
  // set the coupling
  norm(UnitRemoval::E * kappa_ * couplast_);
}

Complex ADDModelFFGGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
					Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}

