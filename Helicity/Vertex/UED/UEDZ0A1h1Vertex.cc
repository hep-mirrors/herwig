// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDZ0A1h1Vertex class.
//

#include "UEDZ0A1h1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig::Helicity;

void UEDZ0A1h1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase << theSin2ThetaW << theKappa;
}

void UEDZ0A1h1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase >> theSin2ThetaW >> theKappa;
  theq2Last = 0.;
  theCoupLast = 0.;
}

ClassDescription<UEDZ0A1h1Vertex> UEDZ0A1h1Vertex::initUEDZ0A1h1Vertex;
// Definition of the static class description member.

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
    if(q2 != theq2Last) {
      theq2Last = q2;
      theCoupLast = 
	theKappa*sqrt(4.*Constants::pi*theUEDBase->alphaEM(q2))/theSin2ThetaW;
    }
    setNorm(theCoupLast); 
  }
  else
    throw HelicityLogicalError() << "UEDZ0A1h1Vertex::setCoupling - "
				 << "There is an unknown particle in this "
				 << "vertex. " << scaA << " " << scaB
				 << Exception::warning;
}
