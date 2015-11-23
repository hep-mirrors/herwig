// -*- C++ -*-
//
// UEDG1G1G0Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDG1G1G0Vertex class.
//

#include "UEDG1G1G0Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDG1G1G0Vertex::UEDG1G1G0Vertex() 
  : theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(1);
  orderInGem(0);
}

void UEDG1G1G0Vertex::doinit() {
  long kkg1 = 5100021;
  addToList(kkg1, kkg1, 21);
  VVVVertex::doinit();
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeNoPIOClass<UEDG1G1G0Vertex,Helicity::VVVVertex>
describeUEDG1G1G0Vertex("Herwig::UEDG1G1G0Vertex", "HwUED.so");

void UEDG1G1G0Vertex::Init() {

  static ClassDocumentation<UEDG1G1G0Vertex> documentation
    ("The UEDG1G1G0Vertex class implements the coupling of the "
     "gluon to two KK excitations of the gluon in the UED model.");

}

void UEDG1G1G0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2, 
				  tcPDPtr part3) {
  long id1(part1->id()), id2(part2->id()), id3(part3->id());
  if( (id1 == ParticleID::g && id2 == 5100021 && id3 == 5100021) ||
      (id2 == ParticleID::g && id1 == 5100021 && id3 == 5100021) ||
      (id3 == ParticleID::g && id1 == 5100021 && id2 == 5100021) ) {
    if(q2 != theq2Last || theCoupLast == 0.) {
      theq2Last = q2;
      theCoupLast = strongCoupling(q2);
    }
    norm(theCoupLast);
  }
  else throw HelicityLogicalError() 
    << "UEDG1G1G0Vertex::setCoupling - "
    << "There is an unknown particle in this vertex "
    << id1 << " " << id2 << " " << id3 << Exception::runerror;
}
