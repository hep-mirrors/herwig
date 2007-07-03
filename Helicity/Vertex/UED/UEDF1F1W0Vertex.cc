// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F1W0Vertex class.
//

#include "UEDF1F1W0Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig::Helicity;

UEDF1F1W0Vertex::UEDF1F1W0Vertex(): theSinThW(0.), theRadius(),
				    theQ2Last(), theCoupLast(0.), 
				    theLeftLast(0.), thefermALast(0),
				    thefermBLast(0) {
  vector<int> anti, ferm, wboson;
  //outgoing W+
  for(unsigned int i = 2; i < 7; i += 2) {
    for(unsigned int j = 1; j < 6; j += 2) {
      anti.push_back(-(5100000 + i));
      ferm.push_back(5100000 + j);
      wboson.push_back(24);
    }
  }
  for(unsigned int i = 2; i < 7; i += 2) {
    for(unsigned int j = 1; j < 6; j += 2) {
      anti.push_back(-(6100000 + i));
      ferm.push_back(6100000 + j);
      wboson.push_back(24);
    }
  }
  for(unsigned int i = 11; i < 17; i += 2) {
    anti.push_back(-(5100001 + i));
    ferm.push_back(5100000 + i);
    wboson.push_back(24);
  }
  //outgoing W-
  for(unsigned int i = 1; i < 6; i += 2) {
     for(unsigned int j = 2 ; j < 7; j += 2) {
       anti.push_back(-(5100000 + i));
       ferm.push_back(5100000 + j);
       wboson.push_back(-24);
     }
  }
  for(unsigned int i = 11; i < 17; i += 2) {
    anti.push_back(-(5100000 + i));
    ferm.push_back(5100001 + i);
    wboson.push_back(-24);
  }
  for(unsigned int i = 1; i < 6; i += 2) {
     for(unsigned int j = 2 ; j < 7; j += 2) {
       anti.push_back(-(6100000 + i));
       ferm.push_back(6100000 + j);
       wboson.push_back(-24);
     }
  }
  //NO MIXING YET
  setList(anti, ferm, wboson);
}


void UEDF1F1W0Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase << theSinThW << ounit(theRadius,1/GeV);
}

void UEDF1F1W0Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase >> theSinThW >> iunit(theRadius,1/GeV);
  theQ2Last = 0.0*GeV2;
  theCoupLast = 0.0;
  theLeftLast = 0.0;
  thefermALast = 0;
  thefermBLast = 0;
}

ClassDescription<UEDF1F1W0Vertex> UEDF1F1W0Vertex::initUEDF1F1W0Vertex;
// Definition of the static class description member.

void UEDF1F1W0Vertex::Init() {

  static ClassDocumentation<UEDF1F1W0Vertex> documentation
    ("This class implements the coupling of a pair of level-1 KK fermions"
     "to an SM W boson");

}

void UEDF1F1W0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long iferm(0), ianti(0);
  if(abs(part1->id()) == ParticleID::Wplus) {
    iferm = part2->id();
    ianti = part3->id();
    if(iferm < 0) swap(iferm, ianti);
  }
  else if(abs(part2->id()) == ParticleID::Wplus) {
    iferm = part1->id();
    ianti = part3->id();
    if(iferm < 0) swap(iferm, ianti);
  }
  else if(abs(part3->id()) == ParticleID::Wplus) {
    iferm = part2->id();
    ianti = part1->id();
    if(iferm < 0) swap(iferm, ianti);
  }
  else
    throw HelicityConsistencyError() << "UEDFFW0Vertex::setCoupling - "
				     << "There is no W boson in this vertex!"
				     << Exception::runerror;
  ianti = abs(ianti);
  bool ferma = (iferm >= 5100001 && iferm <= 5100006) ||
    (iferm >= 6100001 && iferm <= 6100006) || 
    (iferm >= 5100011 && iferm <= 5100016) ||
    (iferm >= 6100011 && iferm <= 6100016); 
  bool fermb = (ianti >= 5100001 && ianti <= 5100006) ||
    (ianti >= 6100001 && ianti <= 6100006) || 
    (ianti >= 5100011 && ianti <= 5100016) ||
    (ianti >= 6100011 && ianti <= 6100016);
  if(ferma && fermb) {
    if(q2 != theQ2Last) {
      theQ2Last = q2;
	theCoupLast = 
	  sqrt(4.*Constants::pi*theUEDBase->alphaEM(q2)/2.)/theSinThW;
    }
    if(iferm != thefermALast || ianti != thefermBLast) {
      thefermALast = iferm;
      thefermBLast = ianti;
      int stateA = ianti/1000000;
      int stateB = iferm/1000000;
      long smID = (stateA == 6) ? iferm - 6100000 : iferm - 5100000;
	// L/R mixing
      double beta = getParticleData(smID)->mass()*theRadius;
      double gamma = beta*beta/(1. + beta*beta);
      double sin2alU = 0.5 - 0.5*sqrt(1. - gamma);
      double cos2alU = 1. - sin2alU; 

      smID = (stateA == 6) ? ianti - 6100000 : ianti - 5100000;
      beta = getParticleData(smID)->mass()*theRadius;
      gamma = beta*beta/(1. + beta*beta);
      double sin2alD = 0.5 - 0.5*sqrt(1. - gamma);
      double cos2alD = 1. - sin2alD;
      if(stateA == 5 && stateB == 5)
	theLeftLast = sqrt(cos2alU*cos2alD);
      else if(stateA == 6 && stateB == 6)
	theLeftLast = sqrt(sin2alU*sin2alD);
      else {}//CHECK MIXING FEYNMAN RULE	
    }
    setNorm(theCoupLast);
    setLeft(theLeftLast);
    setRight(theLeftLast);
  }
  else {
    throw HelicityLogicalError() << "UEDF1F1W0Vertex::setCoupling - "
				 << "There is an unknown particle(s) in the "
				 << "UED F^(1) F^(1) W^(0) vertex. ID: " 
				 << ianti << " " << iferm 
				 << Exception::warning;      
    setNorm(0.0);
    setLeft(0.0);
    setRight(0.0);  
}
}
