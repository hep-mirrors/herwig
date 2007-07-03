// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0Z1Vertex class.
//

#include "UEDF1F0Z1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig::Helicity;

UEDF1F0Z1Vertex::UEDF1F0Z1Vertex() : theSinThetaW(0.), theCosThetaW(0.),
				     theSinThetaOne(0.),
				     theCosWmOne(0.), theq2Last(), 
				     theCoupLast(0.), theKKLast(0), 
				     theLeftLast(0.), theRightLast(0.) {
vector<int> anti, ferm, kkz(42, 5100023);
  //QQ
  for(int i = 1; i < 7; ++i) {
    anti.push_back(-i);
    ferm.push_back(i + 5100000);
    anti.push_back(-(i + 5100000));
    ferm.push_back(i);
    anti.push_back(-i);
    ferm.push_back(i + 6100000);
    anti.push_back(-(i + 6100000));
    ferm.push_back(i);
  }
  //LL
  for(int i = 11; i < 17; ++i) {
    anti.push_back(-i);
    ferm.push_back(i + 5100000);
    anti.push_back(-(i + 5100000));
    ferm.push_back(i);
  }
  for(int i = 11; i < 17; i +=2) {
    anti.push_back(-i);
    ferm.push_back(i + 6100000);
    anti.push_back(-(i + 6100000));
    ferm.push_back(i);
  }
  setList(anti, ferm, kkz);
}

void UEDF1F0Z1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase << theSinThetaW  << theCosThetaW << theSinThetaOne
     << theCosWmOne;
}

void UEDF1F0Z1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase >> theSinThetaW >> theCosThetaW  >> theSinThetaOne
     >> theCosWmOne;
  theq2Last = 0.*GeV2;
  theCoupLast = 0.;
  theKKLast = 0;
  theLeftLast = 0.;
  theRightLast = 0.;
}

ClassDescription<UEDF1F0Z1Vertex> UEDF1F0Z1Vertex::initUEDF1F0Z1Vertex;
// Definition of the static class description member.

void UEDF1F0Z1Vertex::Init() {

  static ClassDocumentation<UEDF1F0Z1Vertex> documentation
    ("This is the coupling of a level-1 KK fermion to an SM fermion and"
     " a level-1 KK Z boson.");

}

void UEDF1F0Z1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long kkparticle(0);
  if(part1->id() == 5100023) {
    kkparticle = (abs(part2->id()) > 5000000) ? abs(part2->id()) : abs(part3->id());
  }
  else if(part2->id() == 5100023) {
    kkparticle = (abs(part1->id()) > 5000000) ? abs(part1->id()) : abs(part3->id());
  }
  else if(part3->id() == 5100023) {
    kkparticle = (abs(part1->id()) > 5000000) ? abs(part1->id()) : abs(part2->id());
  }
  else {
    throw HelicityLogicalError() << "UEDF1F0Z1Vertex::setCoupling - "
				 << "There is no KK Z in this vertex!"
				 << Exception::warning;
    return;
  }
  if( (kkparticle >= 5100001 && kkparticle <= 5100006) ||
      (kkparticle >= 6100001 && kkparticle <= 6100006) ||
      (kkparticle >= 5100011 && kkparticle <= 5100016) ||
      (kkparticle >= 6100011 && kkparticle <= 6100016) ) {
    if(q2 != theq2Last) {
      theq2Last = q2;
      theCoupLast = 
	sqrt(4.*Constants::pi*(theUEDBase->alphaEM(q2)))/theCosThetaW;
    }
    setNorm(theCoupLast);
    if(kkparticle != theKKLast) {
      theKKLast = kkparticle;
      long smID = (kkparticle > 6000000) ? kkparticle - 6100000 
	: kkparticle - 5100000;
      Charge Qf = getParticleData(smID)->charge();
      if(kkparticle/1000000 == 5) {
	double I3f = (abs(smID) % 2 == 0) ? 0.5 : -0.5;
	theLeftLast = (Qf/eplus*theSinThetaOne)
 	  - (I3f*theCosWmOne/theSinThetaW);
	theRightLast = 0.;
      }
      else {
	theLeftLast = 0.;
	theRightLast = Qf/eplus*theSinThetaOne;
      }	
    }
    setLeft(theLeftLast);
    setRight(theRightLast);
  }
  else
    throw HelicityLogicalError() << "UEDF1F0Z1Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex! "
				 << kkparticle
				 << Exception::warning;
} 

