// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDZ0H1H1Vertex class.
//

#include "UEDZ0H1H1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig::Helicity;

void UEDZ0H1H1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase << theSinThetaW << theCosThetaW << theCosTheta2W
     << theMw2 << theR2;
}

void UEDZ0H1H1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase >> theSinThetaW >> theCosThetaW >> theCosTheta2W
     >> theMw2 >> theR2;
  theq2Last = 0.;
  theCoupLast = 0.;
}

ClassDescription<UEDZ0H1H1Vertex> UEDZ0H1H1Vertex::initUEDZ0H1H1Vertex;
// Definition of the static class description member.

void UEDZ0H1H1Vertex::Init() {

  static ClassDocumentation<UEDZ0H1H1Vertex> documentation
    ("This is the coupling of the SM Z boson to the level-1 charged higgs.");

}

void UEDZ0H1H1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long kkhiggs(0);
  if(part1->id() == ParticleID::Z0)
    kkhiggs = abs(part2->id());
  else if(part2->id() == ParticleID::Z0)
    kkhiggs = abs(part1->id());
  else if(part3->id() == ParticleID::Z0)
    kkhiggs = abs(part1->id());
  else {
    throw HelicityLogicalError() << "UEDZ0H1H1Vertex::setCoupling - There is no "
				 << "SM photon in this vertex!." 
				 << Exception::warning;
    return;
  }
  if(kkhiggs == 5100037) {
    if(q2 != theq2Last) {
      theq2Last = q2;
      theCoupLast = 
	Complex(0., 1.)*sqrt(4.*Constants::pi*theUEDBase->alphaEM(q2))/theSinThetaW;
      theCoupLast *= ( (theCosTheta2W/2./theCosThetaW) 
		       - sqr(theCosThetaW)*theMw2*theR2 )/(1. + theMw2*theR2);
    }
    setNorm(theCoupLast);
  }
  else
    throw HelicityLogicalError() << "UEDZ0H1H1Vertex::setCoupling - There is no "
				 << "level-1 higgs in this vertex! " << kkhiggs
				 << Exception::warning;
}
