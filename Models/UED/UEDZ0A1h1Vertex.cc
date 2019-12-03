// -*- C++ -*-
//
// UEDZ0A1h1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDZ0A1h1Vertex class.
//

#include "UEDZ0A1h1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDZ0A1h1Vertex::UEDZ0A1h1Vertex() : theSin2ThetaW(0.), theKappa(0.),	    
				     theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void UEDZ0A1h1Vertex::doinit() {
  addToList(23, 5100036, 5100025);
  VSSVertex::doinit();
  tUEDBasePtr UEDBase = 
    dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!UEDBase)
    throw InitException() << "UEDZ0A1h1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  double sw2 = sin2ThetaW();
  theSin2ThetaW = 2.*sqrt(sw2*(1. - sw2));
  Energy2 mz2 = sqr(getParticleData(23)->mass());
  InvEnergy2 rad2 = sqr(UEDBase->compactRadius());
  theKappa = 1./sqrt(1. + mz2*rad2);
}

void UEDZ0A1h1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theSin2ThetaW << theKappa;
}

void UEDZ0A1h1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theSin2ThetaW >> theKappa;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDZ0A1h1Vertex,VSSVertex>
describeHerwigUEDZ0A1h1Vertex("Herwig::UEDZ0A1h1Vertex", "HwUED.so");

void UEDZ0A1h1Vertex::Init() {

  static ClassDocumentation<UEDZ0A1h1Vertex> documentation
    ("The coupling of an SM Z boson to a level-1 CP-Odd pseudo-scalar "
     "and level 1 higgs.");

}

void UEDZ0A1h1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2, 
				  tcPDPtr part3) {
  long scaA(0), scaB(0);
  if(part1->id() == ParticleID::Z0) {
    scaA = part2->id();
    scaB = part3->id();
  }
  else if(part2->id() == ParticleID::Z0) {
    scaA = part1->id();
    scaB = part3->id();
  }
  else if(part3->id() == ParticleID::Z0) {
    scaA =  part1->id();
    scaB = part2->id();
  }
  else {
    throw HelicityLogicalError() << "UEDZ0A1h1Vertex::setCoupling - "
				 << "There is no SM Z boson in this vertex"
				 << Exception::warning;
  }
  if( (scaA == 5100036 && scaB == 5100025) ||
      (scaB == 5100036 && scaA == 5100025) ) {
    if(q2 != theq2Last || theCoupLast == 0.) {
      theq2Last = q2;
      theCoupLast = theKappa*electroMagneticCoupling(q2)/theSin2ThetaW;
    }
    norm(theCoupLast); 
  }
  else
    throw HelicityLogicalError() << "UEDZ0A1h1Vertex::setCoupling - "
				 << "There is an unknown particle in this "
				 << "vertex. " << scaA << " " << scaB
				 << Exception::warning;
}
