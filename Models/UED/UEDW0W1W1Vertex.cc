// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDW0W1W1Vertex class.
//

#include "UEDW0W1W1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;
UEDW0W1W1Vertex::UEDW0W1W1Vertex() : theZfact(0.), theSinThetaOne(0.), 
				     theq2Last(), theCoupLast(0.) {
  vector<int> first(4), second(4), third(4);
  first[0] = 24;
  second[0] = -5100024;
  third[0] = 5100023;

  first[1] = 5100024;
  second[1] = -24;
  third[1] = 5100023;

  first[2] = 24;
  second[2] = -5100024;
  third[2] = 5100022;

  first[3] = 5100024;
  second[3] = -24;
  third[3] = 5100022;
  setList(first, second, third);
}

void UEDW0W1W1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase << theZfact << theSinThetaOne;
}

void UEDW0W1W1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase  >> theZfact >> theSinThetaOne;
  theq2Last = 0.*GeV2;
  theCoupLast = 0.;
}

ClassDescription<UEDW0W1W1Vertex> UEDW0W1W1Vertex::initUEDW0W1W1Vertex;
// Definition of the static class description member.

void UEDW0W1W1Vertex::Init() {

  static ClassDocumentation<UEDW0W1W1Vertex> documentation
    ("The coupling of an SM W boson to a level 1 KK W and KK Z and KK photon");

}

void UEDW0W1W1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long kkparticle(0);
  double prefact(1.);
  if(abs(part1->id()) == 24) {
    if(part1->id() > 0) prefact = -1.;
    kkparticle = ( abs(part2->id()) == 5100024 ) ? part3->id() : part2->id();
  }
  else if(abs(part2->id()) == 24) {
    if(part2->id() > 0) prefact = -1.;
    kkparticle = ( abs(part1->id()) == 5100024 ) ? part3->id() : part1->id();
  }
  else if(abs(part3->id()) == 24) {
    if(part3->id() > 0) prefact = -1.;
    kkparticle = ( abs(part1->id()) == 5100024 ) ? part2->id() : part1->id();
  }
  else {
    throw HelicityLogicalError() << "UEDW0W1W1Vertex::setCoupling() - "
				 << "There is no W boson in this vertex!"
				 << Exception::warning;
    return;
  }

  if(q2 != theq2Last) {
    theq2Last = q2;
    theCoupLast = sqrt(4.*Constants::pi*theUEDBase->alphaEM(q2));
  }
  if(kkparticle == 5100023)
    setNorm(prefact*theCoupLast*theZfact);
  else if(kkparticle == 5100022)
    setNorm(prefact*theCoupLast*theSinThetaOne);
  else {
    throw HelicityLogicalError() << "UEDW0W1W1Vertex::setCoupling() - "
				 << "There is an unknown particle in this "
				 << "vertex! " << kkparticle
				 << Exception::warning;
    setNorm(0.);
  }
}
