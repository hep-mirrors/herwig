// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDP0H1H1Vertex class.
//

#include "UEDP0H1H1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig::Helicity;

void UEDP0H1H1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase;
}

void UEDP0H1H1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase;
  theq2Last = 0.*GeV2;
  theCoupLast = 0.;
}

ClassDescription<UEDP0H1H1Vertex> UEDP0H1H1Vertex::initUEDP0H1H1Vertex;
// Definition of the static class description member.

void UEDP0H1H1Vertex::Init() {

  static ClassDocumentation<UEDP0H1H1Vertex> documentation
    ("This is the coupling of the SM photon to the level-1 charged higgs.");

}

void UEDP0H1H1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long kkhiggs(0);
  if(part1->id() == ParticleID::gamma)
    kkhiggs = abs(part2->id());
  else if(part2->id() == ParticleID::gamma)
    kkhiggs = abs(part1->id());
  else if(part3->id() == ParticleID::gamma)
    kkhiggs = abs(part1->id());
  else {
    throw HelicityLogicalError() << "UEDP0H1H1Vertex::setCoupling - There is no "
				 << "SM photon in this vertex!." 
				 << Exception::warning;
    return;
  }
  if(kkhiggs == 5100037) {
    if(q2 != theq2Last) {
      theq2Last = q2;
      theCoupLast = 
	Complex(0., 1.)*sqrt(4.*Constants::pi*theUEDBase->alphaEM(q2));
    }
    setNorm(theCoupLast);
  }
  else
    throw HelicityLogicalError() << "UEDP0H1H1Vertex::setCoupling - There is no "
				 << "level-1 higgs in this vertex! " << kkhiggs
				 << Exception::warning;
}
