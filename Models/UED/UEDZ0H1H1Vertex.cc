// -*- C++ -*-
//
// UEDZ0H1H1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDZ0H1H1Vertex class.
//

#include "UEDZ0H1H1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDZ0H1H1Vertex::UEDZ0H1H1Vertex() : theCosThetaW(0.), theCosTheta2W(0.), theMw2(), 
				     theR2(), theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(0);
  orderInGem(1);
}

void UEDZ0H1H1Vertex::doinit() {
  addToList(23, 5100037, -5100037);
  VSSVertex::doinit();
  tUEDBasePtr theUEDBase = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!theUEDBase)
    throw InitException() << "UEDZ0H1H1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theCosThetaW = sqrt(1. - sin2ThetaW());
  theCosTheta2W = 1. - 2.*sin2ThetaW();
  theMw2 = sqr(getParticleData(24)->mass());
  theR2 = sqr(theUEDBase->compactRadius());
}

void UEDZ0H1H1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theCosThetaW << theCosTheta2W
     << ounit(theMw2,GeV2) << ounit(theR2,1/GeV2);
}

void UEDZ0H1H1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theCosThetaW >> theCosTheta2W
     >> iunit(theMw2,GeV2) >> iunit(theR2,1/GeV2);
}

ClassDescription<UEDZ0H1H1Vertex> UEDZ0H1H1Vertex::initUEDZ0H1H1Vertex;
// Definition of the static class description member.

void UEDZ0H1H1Vertex::Init() {

  static ClassDocumentation<UEDZ0H1H1Vertex> documentation
    ("This is the coupling of the SM Z boson to the level-1 charged higgs.");

}

void UEDZ0H1H1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long kkhiggs(0);
  if(part1->id() == ParticleID::Z0)
    kkhiggs = abs(part2->id());
  else if(part2->id() == ParticleID::Z0)
    kkhiggs = abs(part1->id());
  else if(part3->id() == ParticleID::Z0)
    kkhiggs = abs(part1->id());
  else {
    throw HelicityLogicalError() << "UEDZ0H1H1Vertex::setCoupling - There is no "
				 << "SM photon in this vertex!." 
				 << Exception::warning;
    return;
  }
  if(kkhiggs == 5100037) {
    if(q2 != theq2Last || theCoupLast == 0.) {
      theq2Last = q2;
      theCoupLast = 
	Complex(0., 1.)*weakCoupling(q2);
      theCoupLast *= ( (theCosTheta2W/2./theCosThetaW) 
		       - sqr(theCosThetaW)*theMw2*theR2 )/(1. + theMw2*theR2);
    }
    norm(theCoupLast);
  }
  else
    throw HelicityLogicalError() << "UEDZ0H1H1Vertex::setCoupling - There is no "
				 << "level-1 higgs in this vertex! " << kkhiggs
				 << Exception::warning;
}
