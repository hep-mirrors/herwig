// -*- C++ -*-
//
// UEDF1F0G1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0G1Vertex class.
//

#include "UEDF1F0G1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F0G1Vertex::UEDF1F0G1Vertex() : theq2Last(ZERO), theCoupLast(0.) {
  orderInGs(1);
  orderInGem(0);
}

NoPIOClassDescription<UEDF1F0G1Vertex> UEDF1F0G1Vertex::initUEDF1F0G1Vertex;
// Definition of the static class description member.

void UEDF1F0G1Vertex::Init() {

  static ClassDocumentation<UEDF1F0G1Vertex> documentation
    ("This class implements the F^1-F^0-G^1 vertex.");

}

void UEDF1F0G1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long ifermN;
  if(part1->id() == 5100021) {
    if(abs(part2->id()) > 5100000)
      ifermN = part2->id();
    else
      ifermN = part3->id();
  }
  else if(part2->id() == 5100021) {
    if(abs(part1->id()) > 5100000)
      ifermN = part1->id();
    else
      ifermN = part3->id();
  }
  else if(part3->id() == 5100021) {
    if(abs(part2->id()) > 5100000)
      ifermN = part2->id();
    else
      ifermN = part1->id();
  }
  else
    throw HelicityLogicalError() << "UEDF1F0G1Vertex::setCoupling - "
				 << "There is no KK gluon in this vertex!"
				 << Exception::warning;
  if((abs(ifermN) >= 5100001 && abs(ifermN) <= 5100006) ||
     (abs(ifermN) >= 6100001 && abs(ifermN) <= 6100006)) {
    if(q2 != theq2Last || theCoupLast == 0.) {
      theq2Last = q2;
      theCoupLast = -strongCoupling(q2);
    }
    norm(theCoupLast);
    int state = abs(ifermN)/1000000;
    if(state == 5) {
      left(1.);
      right(0.);
    }
    else {
      left(0.);
      right(1.);
    }
  }
  else
    throw HelicityLogicalError() << "UEDF1F0G1Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex! "
				 << ifermN
				 << Exception::warning;
}

void UEDF1F0G1Vertex::doinit() {
  long boson = 5100021;
  //QQ
  for(long i = 1; i < 7; ++i) {
    addToList(-i, i + 5100000, boson);
    addToList(-(i + 5100000), i, boson);

    addToList(-i, i + 6100000, boson);
    addToList(-(i + 6100000), i, boson);
  }
  FFVVertex::doinit();
}
