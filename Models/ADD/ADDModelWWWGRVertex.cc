// -*- C++ -*-
//
// ADDModelWWWGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelWWWGRVertex class.
//

#include "ADDModelWWWGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

ADDModelWWWGRVertex::ADDModelWWWGRVertex() 
  : kappa_(ZERO), r_(ZERO), couplast_(0.),
    q2last_(ZERO), zfact_(0.) {
  // order in the couplings
  orderInGem(2);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

void ADDModelWWWGRVertex::doinit() {
  addToList(24,-24, 22, 39);
  addToList(24,-24, 23, 39);
  VVVTVertex::doinit();
  zfact_ = sqrt((1.-sin2ThetaW())/sin2ThetaW());
  // set the graviton coupling 
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD) 
    throw Exception() 
      << "Must have ADDModel in ADDModelWWWGRVertex::doinit()"
      << Exception::runerror;
  kappa_=2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

void ADDModelWWWGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << zfact_ << ounit(r_,GeV);
}
void ADDModelWWWGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> zfact_ >> iunit(r_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModelWWWGRVertex,VVVTVertex>
describeHerwigADDModelWWWGRVertex("Herwig::ADDModelWWWGRVertex", "HwADDModel.so");

void ADDModelWWWGRVertex::Init() {
 static ClassDocumentation<ADDModelWWWGRVertex> documentation
    ("The ADDModelWWWGRVertex class is the four point coupling"
     " of three vector bosons and a graviton in the Randell-Sundrum model.");
  
}


// couplings for the WWWGR vertex
void ADDModelWWWGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
				     tcPDPtr c, tcPDPtr) {
  int ida=a->id();
  int idb=b->id();
  int idc=c->id();
  // first the overall normalisation
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = electroMagneticCoupling(q2);
    q2last_=q2;
  }
  // W- W+ photon and cylic perms
  if((ida==-24 && idb== 24 && idc== 22) ||
     (ida== 22 && idb==-24 && idc== 24) || 
     (ida== 24 && idb== 22 && idc==-24) )
    norm(Complex(couplast_*kappa_*UnitRemoval::E));
  // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) ||
	  (ida== 22 && idb== 24 && idc==-24) || 
	  (ida==-24 && idb== 22 && idc== 24) )
    norm(-Complex(couplast_*kappa_*UnitRemoval::E));
  // W- W+ Z and cylic perms
  else if((ida==-24 && idb== 24 && idc== 23) ||
	  (ida== 23 && idb==-24 && idc== 24) || 
	  (ida== 24 && idb== 23 && idc==-24) )
    norm(Complex(couplast_*zfact_*kappa_*UnitRemoval::E));
  // W+ W- Z (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 23) ||
	  (ida== 23 && idb== 24 && idc==-24) || 
	  (ida==-24 && idb== 23 && idc== 24) )
    norm(-Complex(couplast_*zfact_*kappa_*UnitRemoval::E));
  else assert(false);
}

Complex ADDModelWWWGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
					Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}
