// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDG1G1G0Vertex class.
//

#include "UEDG1G1G0Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig::Helicity;

void UEDG1G1G0Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase;
}

void UEDG1G1G0Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase;
  theq2Last = 0.;
  theCoupLast = 0.;
}

ClassDescription<UEDG1G1G0Vertex> UEDG1G1G0Vertex::initUEDG1G1G0Vertex;
// Definition of the static class description member.

void UEDG1G1G0Vertex::Init() {

  static ClassDocumentation<UEDG1G1G0Vertex> documentation
    ("There is no documentation for the UEDG1G1G0Vertex class");

}

void UEDG1G1G0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2, 
				  tcPDPtr part3) {
  long id1(part1->id()), id2(part2->id()), id3(part3->id());
  if( (id1 == ParticleID::g && id2 == 5100021 && id3 == 5100021) ||
      (id2 == ParticleID::g && id1 == 5100021 && id3 == 5100021) ||
      (id3 == ParticleID::g && id1 == 5100021 && id2 == 5100021) ) {
    if(q2 != theq2Last) {
      theq2Last = q2;
      theCoupLast = sqrt(4.*Constants::pi*theUEDBase->alphaS(q2));
    }
    setNorm(theCoupLast);
  }
  else {
    throw HelicityLogicalError() << "UEDG1G1G0Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex "
				 << id1 << " " << id2 << " " << id3 
				 << Exception::warning;
    setNorm(0.);
  }
}
