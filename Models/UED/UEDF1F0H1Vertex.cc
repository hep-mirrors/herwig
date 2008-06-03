// -*- C++ -*-
//
// UEDF1F0H1Vertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0H1Vertex class.
//

#include "UEDF1F0H1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F0H1Vertex::UEDF1F0H1Vertex() : theRadius(), theMw(), 
				     theSinThetaW(0.), theq2Last(-1.*GeV2),
				     theCoupLast(0.), theLeftLast(0.),
				     theRightLast(0.), theKKLast(0),
				     theSMLast(0) {
  vector<long> anti, ferm, kkhiggs;
  for(long i = 1; i < 6; i += 2) {
    //outgoing H+
    anti.push_back(-(5100001 + i));
    ferm.push_back(i);
    kkhiggs.push_back(5100037);

    anti.push_back(-(6100001 + i));
    ferm.push_back(i);
    kkhiggs.push_back(5100037);

    anti.push_back(-i - 1);
    ferm.push_back(5100000 + i);
    kkhiggs.push_back(5100037);

    anti.push_back(-i - 1);
    ferm.push_back(6100000 + i);
    kkhiggs.push_back(5100037);
    
    //outgoing H-
    anti.push_back(-5100000 - i);
    ferm.push_back(i + 1);
    kkhiggs.push_back(-5100037);
    
    anti.push_back(-6100000 - i);
    ferm.push_back(i + 1);
    kkhiggs.push_back(-5100037);    
    
    anti.push_back(-i);
    ferm.push_back(5100001 + i);
    kkhiggs.push_back(-5100037);
  }
  setList(anti, ferm, kkhiggs);
}

void UEDF1F0H1Vertex::doinit() throw(InitException) {
  FFSVertex::doinit();
  tUEDBasePtr UEDBase = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!UEDBase)
    throw InitException() << "UEDF1F0H1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theRadius = UEDBase->compactRadius();
  theSinThetaW = sqrt(UEDBase->sin2ThetaW());
  theMw = getParticleData(24)->mass();
  orderInGs(0);
  orderInGem(1);
}


void UEDF1F0H1Vertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theRadius,1/GeV) << ounit(theMw,GeV) << theSinThetaW ;
}

void UEDF1F0H1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theRadius,1/GeV) >> iunit(theMw,GeV) >> theSinThetaW;
}

ClassDescription<UEDF1F0H1Vertex> UEDF1F0H1Vertex::initUEDF1F0H1Vertex;
// Definition of the static class description member.

void UEDF1F0H1Vertex::Init() {

  static ClassDocumentation<UEDF1F0H1Vertex> documentation
    ("The coupling of a KK-1 fermion to an SM fermion and a charged higgs.");

}

void UEDF1F0H1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3, int) {
  long kkferm(0), smferm(0), kkhiggs(0);
  if(abs(part1->id()) == 5100037) {
    kkhiggs = part1->id();
    if(abs(part2->id()) > 5100000) {
      kkferm = abs(part2->id());
      smferm = part3->id();
    }
    else {
      kkferm = abs(part3->id());
      smferm = part2->id();
    }
  }
  else if(abs(part2->id()) == 5100037) {
    kkhiggs = part2->id();
    if(abs(part1->id()) > 5100000) {
      kkferm = abs(part1->id());
      smferm = part3->id();
    }
    else {
      kkferm = abs(part3->id());
      smferm = part1->id();
    }
  }
  else if(abs(part3->id()) == 5100037) {
    kkhiggs = part3->id();
    if(abs(part1->id()) > 5100000) {
      kkferm = abs(part1->id());
      smferm = part2->id();
    }
    else {
      kkferm = abs(part2->id());
      smferm = part1->id();
    }
  }
  else {
    throw HelicityLogicalError() << "UEDF1F0H1Vertex::setCoupling - There is "
				 << "no KK1 charged Higgs in this vertex."
				 << Exception::warning;
    return;
  }
  bool cc;
  if(kkhiggs > 0)
    cc = (smferm > 0) ? false : true;
  else
    cc = (smferm > 0) ? true : false;
  smferm = abs(smferm);
  if( (smferm >= 1 && smferm <= 6) ||
      (smferm >= 11 && smferm <= 16) ) {
    if(q2 != theq2Last) {
      theq2Last = q2;
      theCoupLast = sqrt(0.5)*weakCoupling(q2);
      theCoupLast *= theRadius/sqrt(1. + sqr(theMw*theRadius)) * UnitRemoval::E;
    }
    setNorm(theCoupLast);
    if(kkferm != theKKLast || abs(smferm) != theSMLast) {
      long smID = (kkferm > 6100000) ? kkferm - 6100000 : kkferm - 5100000;
      Energy smMass = getParticleData(smID)->mass();
      double beta = smMass*theRadius;
      double gamma = beta*beta/(1. + beta*beta);
      double sinAl = sqrt(0.5 - 0.5*sqrt(1. - gamma));
      double cosAl = sqrt(1. - sinAl*sinAl);
      unsigned int kkstate = kkferm/1000000;

      theRightLast = smMass/theRadius/theMw * UnitRemoval::InvE;
      if(kkstate == 5) {
	theLeftLast = (theMw*cosAl - (smMass*sinAl/theRadius/theMw)) 
	  * UnitRemoval::InvE;
	theRightLast *= cosAl;
      }
      else {
	theLeftLast = (theMw*sinAl + (smMass*cosAl/theRadius/theMw))
	  * UnitRemoval::InvE;
	theRightLast *= -sinAl;
      }
    }
    if(kkhiggs < 0) {
      setLeft(theRightLast);
      setRight(theLeftLast);
    }
    else {
      setLeft(theLeftLast);
      setRight(theRightLast);
    }
  }
  else 
    throw HelicityLogicalError() << "UEDF1F0H1Vertex::setCoupling - There is an "
				 << "unknown particle in this vertex " 
				 << smferm << " " << kkferm
				 << Exception::warning;
  
}
