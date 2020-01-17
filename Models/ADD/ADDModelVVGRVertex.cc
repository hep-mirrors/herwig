// -*- C++ -*-
//
// ADDModelVVGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelVVGRVertex class.
//

#include "ADDModelVVGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void ADDModelVVGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << ounit(r_,GeV);
}

void ADDModelVVGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> iunit(r_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModelVVGRVertex,VVTVertex>
describeHerwigADDModelVVGRVertex("Herwig::ADDModelVVGRVertex", "HwADDModel.so");

void ADDModelVVGRVertex::Init() {
 static ClassDocumentation<ADDModelVVGRVertex> documentation
    ("The ADDModelVVGRVertex class is the implementation"
     " of the ADDModel vector-vector-graviton vertex");
  
}
  
void ADDModelVVGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
  norm(Complex(UnitRemoval::E * kappa_));
}

ADDModelVVGRVertex::ADDModelVVGRVertex() : kappa_(ZERO), r_(ZERO) {
  orderInGem(1);
  orderInGs (0);
  colourStructure(ColourStructure::DELTA);
}

void ADDModelVVGRVertex::doinit() {
  addToList(23,23,39);
  addToList(22,22,39);
  addToList(24,-24,39);
  addToList(21,21,39);
  VVTVertex::doinit();
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD)
    throw Exception() << "Must be ADDModel in ADDModelVVGRVertex::doinit()"
		      << Exception::runerror;
  kappa_=2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

Complex ADDModelVVGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
				       Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}
