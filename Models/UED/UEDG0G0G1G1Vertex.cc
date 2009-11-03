// -*- C++ -*-
//
// UEDG0G0G1G1Vertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDG0G0G1G1Vertex class.
//

#include "UEDG0G0G1G1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

UEDG0G0G1G1Vertex::UEDG0G0G1G1Vertex() : 
  theq2Last(ZERO), theCoupLast(0.) {
  vector<long> kk1g(1, 5100021), smgl(1, 21);
  setList(smgl, smgl, kk1g, kk1g);
}

void UEDG0G0G1G1Vertex::doinit() {
  VVVVVertex::doinit();
  orderInGs(2);
  orderInGem(0);
}

NoPIOClassDescription<UEDG0G0G1G1Vertex> UEDG0G0G1G1Vertex::initUEDG0G0G1G1Vertex;
// Definition of the static class description member.

void UEDG0G0G1G1Vertex::Init() {

  static ClassDocumentation<UEDG0G0G1G1Vertex> documentation
    ("This class implements the coupling of a pair of SM gluons to"
     "a pair of UED level-1 KK gluons.");

}

void UEDG0G0G1G1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2, 
				    tcPDPtr part3, tcPDPtr part4) {
  int ismg(0), ikkg(0);
  vector<tcPDPtr> particles(4);
  particles[0] = part1; particles[1] = part2;
  particles[2] = part3; particles[3] = part4;
  for(vector<long>::size_type i = 0; i < 4; ++i) {
    if(particles[i]->id() == ParticleID::g) ++ismg;
    if(particles[i]->id() == 5100021) ++ikkg;
  }
  if(ismg == 2 && ikkg == 2) { 
    if(q2 != theq2Last || theCoupLast == 0. ) {
      theq2Last = q2;
      theCoupLast = sqr(strongCoupling(q2));
    }
    setNorm(theCoupLast);
    setType(1); setOrder(0,1,2,3);
  }
  else {
    throw HelicityLogicalError() << "UEDG0G0G1G1Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex "
				 << part1->id() << " " << part2->id() << " " 
				 << part3->id() << " " << part4->id()
				 << Exception::warning;
    setNorm(0.);
  }
}
