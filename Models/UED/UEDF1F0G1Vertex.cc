// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0G1Vertex class.
//

#include "UEDF1F0G1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F0G1Vertex::UEDF1F0G1Vertex() {
  vector<int> anti, ferm, boson(24, 5100021);
  //QQ
    for(int i = 1; i < 7; ++i) {
      anti.push_back(-i);
      ferm.push_back(i + 5100000);
      anti.push_back(-(i + 5100000));
      ferm.push_back(i);
    }
    for(int i = 1; i < 7; ++i) {
      anti.push_back(-i);
      ferm.push_back(i + 6100000);
      anti.push_back(-(i + 6100000));
      ferm.push_back(i);
      
    }
    setList(anti, ferm, boson);
}

void UEDF1F0G1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase;
}

void UEDF1F0G1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase;
  theq2Last = 0.*GeV2;
  theCoupLast = 0.;
}

ClassDescription<UEDF1F0G1Vertex> UEDF1F0G1Vertex::initUEDF1F0G1Vertex;
// Definition of the static class description member.

void UEDF1F0G1Vertex::Init() {

  static ClassDocumentation<UEDF1F0G1Vertex> documentation
    ("This class implements the F^1-F^0-G^1 vertex.");

}

void UEDF1F0G1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long ifermN;
  if(part1->id() == 5100021) {
    if(abs(part2->id()) > 5100000)
      ifermN = part2->id();
    else
      ifermN = part3->id();
  }
  else if(part2->id() == 5100021) {
    if(abs(part1->id()) > 5100000)
      ifermN = part1->id();
    else
      ifermN = part3->id();
  }
  else if(part3->id() == 5100021) {
    if(abs(part2->id()) > 5100000)
      ifermN = part2->id();
    else
      ifermN = part1->id();
  }
  else
    throw HelicityLogicalError() << "UEDF1F0G1Vertex::setCoupling - "
				 << "There is no KK gluon in this vertex!"
				 << Exception::warning;
  if((abs(ifermN) >= 5100001 && abs(ifermN) <= 5100006) ||
     (abs(ifermN) >= 6100001 && abs(ifermN) <= 6100006)) {
    if(q2 != theq2Last) {
      theq2Last = q2;
      theCoupLast = -sqrt(4.*Constants::pi*(theUEDBase->alphaS(q2)));
    }
    setNorm(theCoupLast);
    int state = abs(ifermN)/1000000;
    if(state == 5) {
      setLeft(1.);
      setRight(0.);
    }
    else {
      setLeft(0.);
      setRight(1.);
    }
  }
  else
    throw HelicityLogicalError() << "UEDF1F0G1Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex! "
				 << ifermN
				 << Exception::warning;
}
