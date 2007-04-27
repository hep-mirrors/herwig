// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0H1Vertex class.
//

#include "UEDF1F0H1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig::Helicity;

UEDF1F0H1Vertex::UEDF1F0H1Vertex() : theRadius(0.), theMw(0.), 
				     theSinThetaW(0.), theq2Last(0.),
				     theCoupLast(0.), theLeftLast(0.),
				     theRightLast(0.), theKKLast(0),
				     theSMLast(0) {
  vector<int> anti, ferm, kkhiggs;
  for(unsigned int i = 1; i < 6; i += 2) {
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


void UEDF1F0H1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase << theRadius << theMw << theSinThetaW ;
}

void UEDF1F0H1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase >> theRadius >> theMw >> theSinThetaW;
  theq2Last = 0.;
  theCoupLast = 0.;
  theLeftLast = 0.;
  theRightLast = 0.;
  theKKLast = 0;
  theSMLast = 0;
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
      theCoupLast = sqrt(2.*Constants::pi*theUEDBase->alphaEM(q2))/theSinThetaW;
      theCoupLast *= theRadius/sqrt(1. + theMw*theMw*theRadius*theRadius);
    }
    setNorm(theCoupLast);
    if(kkferm != theKKLast || abs(smferm) != theSMLast) {
      long smID = (kkferm > 6100000) ? kkferm - 6100000 : kkferm - 5100000;
      Energy smMass = getParticleData(smID)->mass();
      smMass = 0.;
      double beta = smMass*theRadius;
      double gamma = beta*beta/(1. + beta*beta);
      double sinAl = sqrt(0.5 - 0.5*sqrt(1. - gamma));
      double cosAl = sqrt(1. - sinAl*sinAl);
      unsigned int kkstate = kkferm/1000000;

      theRightLast = smMass/theRadius/theMw;
      if(kkstate == 5) {
	theLeftLast = theMw*cosAl - (smMass*sinAl/theRadius/theMw);
	theRightLast = cosAl;
      }
      else {
	theLeftLast = theMw*sinAl + (smMass*cosAl/theRadius/theMw);
	theRightLast =-sinAl;
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
