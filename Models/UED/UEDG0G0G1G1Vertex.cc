// -*- C++ -*-
//
// UEDG0G0G1G1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDG0G0G1G1Vertex class.
//

#include "UEDG0G0G1G1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "UEDBase.h"

using namespace Herwig;

UEDG0G0G1G1Vertex::UEDG0G0G1G1Vertex() : 
  theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(2);
  orderInGem(0);
}

void UEDG0G0G1G1Vertex::doinit() {
  long kk1g = 5100021, smgl = 21;
  addToList(smgl, smgl, kk1g, kk1g);
  VVVVVertex::doinit();
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeNoPIOClass<UEDG0G0G1G1Vertex,Helicity::VVVVVertex>
describeUEDG0G0G1G1Vertex("Herwig::UEDG0G0G1G1Vertex", "HwUED.so");

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
  assert(ismg == 2 && ikkg == 2);
  if(q2 != theq2Last || theCoupLast == 0. ) {
    theq2Last = q2;
    theCoupLast = sqr(strongCoupling(q2));
  }
  norm(theCoupLast);
  setType(1);
  setOrder(0,1,2,3);
}
