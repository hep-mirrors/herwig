// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0W1Vertex class.
//

#include "UEDF1F0W1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig::Helicity;

UEDF1F0W1Vertex::UEDF1F0W1Vertex() : theSinThetaW2(0.0), theq2Last(0.0),
				     theCoupLast(0.0) {
  vector<int> ferm, anti, wboson;
  //outgoing W+
  for(unsigned int i = 2; i < 7; i += 2) {
    for(unsigned int j = 1; j < 6; j += 2) {
      anti.push_back(-i);
      ferm.push_back(5100000 + j);
      wboson.push_back(5100024);
      anti.push_back(-(5100000 + i));
      ferm.push_back(j);
      wboson.push_back(5100024);
    }
  }
  for(unsigned int i = 11; i < 17; i += 2) {
    anti.push_back(-i - 1);
    ferm.push_back(5100000 + i);
    wboson.push_back(5100024);
    anti.push_back(-(5100001 + i));
    ferm.push_back(i);
    wboson.push_back(5100024);
  }
  //outgoing W-
  for(unsigned int i = 1; i < 6; i += 2) {
    for(unsigned int j = 2 ; j < 7; j += 2) {
      anti.push_back(-i);
      ferm.push_back(5100000 + j);
      wboson.push_back(-5100024);
      anti.push_back(-(5100000 + i));
      ferm.push_back(j);
      wboson.push_back(-5100024);
    }
  }
  for(unsigned int i = 11; i < 17; i += 2) {
    anti.push_back(-i);
    ferm.push_back(5100001 + i);
    wboson.push_back(-5100024);
    anti.push_back(-(5100000 + i));
    ferm.push_back(i + 1);
    wboson.push_back(-5100024);
  }
  setList(anti, ferm, wboson);
}

void UEDF1F0W1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase << theSinThetaW2;
}

void UEDF1F0W1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase >> theSinThetaW2;
}

ClassDescription<UEDF1F0W1Vertex> UEDF1F0W1Vertex::initUEDF1F0W1Vertex;
// Definition of the static class description member.

void UEDF1F0W1Vertex::Init() {

  static ClassDocumentation<UEDF1F0W1Vertex> documentation
    ("This is the coupling of a KK1 W boson to a KK1 fermion and "
     "a SM fermion.");

}

void UEDF1F0W1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long kkparticle(0);
  if(abs(part1->id()) == 5100024) {
    kkparticle = (abs(part2->id()) > 5000000) ? abs(part2->id()) : abs(part3->id());
  }
  else if(abs(part2->id()) == 5100024) {
    kkparticle = (abs(part1->id()) > 5000000) ? abs(part1->id()) : abs(part3->id());
  }
  else if(abs(part3->id()) == 5100024) {
    kkparticle = (abs(part1->id()) > 5000000) ? abs(part1->id()) : abs(part2->id());
  }
  else {
    throw HelicityLogicalError() << "UEDF1F0W1Vertex::setCoupling - "
				 << "There is no KK W in this vertex!"
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
	-sqrt(4.*Constants::pi*theUEDBase->alphaEM(q2)/2./theSinThetaW2);
    }
    setNorm(theCoupLast);
    setLeft(1.);
    setRight(0.);
  }
  else
    throw HelicityLogicalError() << "UEDF1F0W1Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex! "
				 << kkparticle
				 << Exception::warning;
}

