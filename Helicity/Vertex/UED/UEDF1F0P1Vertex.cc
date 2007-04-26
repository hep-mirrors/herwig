// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0P1Vertex class.
//

#include "UEDF1F0P1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig::Helicity;

UEDF1F0P1Vertex::UEDF1F0P1Vertex() : theSinThetaW(0.), theCosThetaW(0.),
				     theCosThetaOne(0.), theSinWmOne(0.),
				     theq2Last(0.), theCoupLast(0.) ,
				     theKKLast(0), theLeftLast(0.), 
				     theRightLast(0.) {
  vector<int> anti, ferm, kkphot(36, 5100022);
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
  for(int i = 11; i < 17; i += 2) {
    anti.push_back(-i);
    ferm.push_back(i + 5100000);
    anti.push_back(-(i + 5100000));
    ferm.push_back(i);
    anti.push_back(-i);
    ferm.push_back(i + 6100000);
    anti.push_back(-(i + 6100000));
    ferm.push_back(i);
  }
  setList(anti, ferm, kkphot);
}


void UEDF1F0P1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase << theSinThetaW << theCosThetaW << theCosThetaOne
     << theSinWmOne;
}

void UEDF1F0P1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase >> theSinThetaW >> theCosThetaW >> theCosThetaOne
     >> theSinWmOne;
  theq2Last = 0.;
  theCoupLast = 0.;
  theKKLast = 0;
  theLeftLast = 0.;
  theRightLast = 0.;
}

ClassDescription<UEDF1F0P1Vertex> UEDF1F0P1Vertex::initUEDF1F0P1Vertex;
// Definition of the static class description member.

void UEDF1F0P1Vertex::Init() {

  static ClassDocumentation<UEDF1F0P1Vertex> documentation
    ("This is the implementation of the level-1 KK fermion, SM fermion "
     "and level-1 KK photon");

}

void UEDF1F0P1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long kkparticle(0);
  if(part1->id() == 5100022) {
    kkparticle = (abs(part2->id()) > 5000000) ? abs(part2->id()) : abs(part3->id());
  }
  else if(part2->id() == 5100022) {
    kkparticle = (abs(part1->id()) > 5000000) ? abs(part1->id()) : abs(part3->id());
  }
  else if(part3->id() == 5100022) {
    kkparticle = (abs(part1->id()) > 5000000) ? abs(part1->id()) : abs(part2->id());
  }
  else {
    throw HelicityLogicalError() << "UEDF1F0P1Vertex::setCoupling - "
				 << "There is no KK photon in this vertex!"
				 << Exception::warning;
    return;
  }
  if( (kkparticle >= 5100001 && kkparticle <= 5100006) ||
      (kkparticle >= 6100001 && kkparticle <= 6100006) ||
      (kkparticle >= 5100011 && kkparticle <= 5100016) ||
      (kkparticle >= 6100011 && kkparticle <= 6100016) ) {
    if(q2 != theq2Last) {
      theq2Last = q2;
      theCoupLast = -sqrt(4.*Constants::pi*(theUEDBase->alphaEM(q2)))/theCosThetaW;
    }
    setNorm(theCoupLast);
    if(kkparticle != theKKLast) {
      theKKLast = kkparticle;
      long smID = (kkparticle > 6000000) ? kkparticle - 6100000 
	: kkparticle - 5100000;
      Charge Qf = getParticleData(smID)->charge();
      if(kkparticle/1000000 == 5) {
	double I3f = (smID % 2 == 0) ? 0.5 : -0.5;
	theLeftLast = Qf*theCosThetaOne - (I3f*theSinWmOne/theSinThetaW);
	theRightLast = 0.;
      }
      else {
	theLeftLast = 0.;
	theRightLast = Qf*theCosThetaOne;
      }
    }
    setLeft(theLeftLast);
    setRight(theRightLast);
  }
  else 
    throw HelicityLogicalError() << "UEDF1F0P1Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex! "
				 << kkparticle
				 << Exception::warning;
}
