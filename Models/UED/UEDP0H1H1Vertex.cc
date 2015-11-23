// -*- C++ -*-
//
// UEDP0H1H1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDP0H1H1Vertex class.
//

#include "UEDP0H1H1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

UEDP0H1H1Vertex::UEDP0H1H1Vertex() : theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(0);
  orderInGem(1);
}

void UEDP0H1H1Vertex::doinit() {
  addToList(22, 5100037, -5100037);
  VSSVertex::doinit();
}

NoPIOClassDescription<UEDP0H1H1Vertex> UEDP0H1H1Vertex::initUEDP0H1H1Vertex;
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
    if(q2 != theq2Last || theCoupLast == 0.) {
      theq2Last = q2;
      theCoupLast = Complex(0., 1.)*electroMagneticCoupling(q2);
    }
    norm(theCoupLast);
  }
  else
    throw HelicityLogicalError() << "UEDP0H1H1Vertex::setCoupling - There is no "
				 << "level-1 higgs in this vertex! " << kkhiggs
				 << Exception::warning;
}
