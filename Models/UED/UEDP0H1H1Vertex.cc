// -*- C++ -*-
//
// UEDP0H1H1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDP0H1H1Vertex class.
//

#include "UEDP0H1H1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

UEDP0H1H1Vertex::UEDP0H1H1Vertex() : theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void UEDP0H1H1Vertex::doinit() {
  addToList(22, 5100037, -5100037);
  VSSVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<UEDP0H1H1Vertex,VSSVertex>
describeHerwigUEDP0H1H1Vertex("Herwig::UEDP0H1H1Vertex", "HwUED.so");

void UEDP0H1H1Vertex::Init() {

  static ClassDocumentation<UEDP0H1H1Vertex> documentation
    ("This is the coupling of the SM photon to the level-1 charged higgs.");

}
#ifndef NDEBUG
void UEDP0H1H1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr ) {
#else
void UEDP0H1H1Vertex::setCoupling(Energy2 q2, tcPDPtr , tcPDPtr part2,
				  tcPDPtr ) {
#endif

  assert(part1->id()==ParticleID::gamma);
  assert(abs(part2->id()) == 5100037);
  if(q2 != theq2Last || theCoupLast == 0.) {
    theq2Last = q2;
    theCoupLast = electroMagneticCoupling(q2);
  }
  if(part2->id()>0) 
    norm(-theCoupLast);
  else
    norm( theCoupLast);
}
