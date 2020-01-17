// -*- C++ -*-
//
// UEDF1F1G0Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F1G0Vertex class.
//

#include "UEDF1F1G0Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F1G0Vertex::UEDF1F1G0Vertex() 
  : theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TFUND);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<UEDF1F1G0Vertex,FFVVertex>
describeHerwigUEDF1F1G0Vertex("Herwig::UEDF1F1G0Vertex", "HwUED.so");

void UEDF1F1G0Vertex::Init() {

  static ClassDocumentation<UEDF1F1G0Vertex> documentation
    ("This class implements the F^1 F^1 G^0 vertex.");

}

void UEDF1F1G0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long iferm;
  if(part1->id() == ParticleID::g)
    iferm = abs(part2->id());
  else if(part2->id() == ParticleID::g)
    iferm = abs(part1->id());
  else if(part3->id() == ParticleID::g)
    iferm = abs(part1->id());
  else
    throw HelicityLogicalError() << "UEDF1F1G0Vertex::setCoupling - "
				 << "There is no gluon in this vertex!"
				 << Exception::warning;
  if((iferm >= 5100001 && iferm <= 5100006) ||
     (iferm >= 6100001 && iferm <= 6100006)) {
    if(q2 != theq2Last || theCoupLast ==0. ) {
      theCoupLast = -strongCoupling(q2);
      theq2Last=q2;
    }
    norm(theCoupLast);
    left(1.);
    right(1.);
  }
  else
    throw HelicityLogicalError() << "UEDF1F1G0Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex! "
				 << iferm
				 << Exception::warning;
}
void UEDF1F1G0Vertex::doinit() {
  long boson = 21;
  //QQ
  for(long i = 5100001; i < 6100007; ++i) {
    if(i == 5100007) i += 999994;
    addToList(-i, i, boson);
  }
  FFVVertex::doinit();
}
