// -*- C++ -*-
//
// UEDF1F1P0Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F1P0Vertex class.
//

#include "UEDF1F1P0Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F1P0Vertex::UEDF1F1P0Vertex() : theCoupLast(0.0), theq2Last(ZERO),
				     thefermLast(0), theLRLast(0.0), 
				     theCharges(3) {
  orderInGs(0);
  orderInGem(1);
}

void UEDF1F1P0Vertex::doinit() {
  long photon = 22;
  //quarks
  for(long i = 1; i < 7; ++i) {
    //left
    addToList(-5100000 - i, 5100000 + i, photon);
    //right
    addToList(-6100000 - i, 6100000 + i, photon);
  }
  //leptons
  for(long i = 11; i < 17; i += 2) {
    //left
    addToList(-5100000 - i, 5100000 + i, photon);
    //right
    addToList(-6100000 - i, 6100000 + i, photon);
  }
  FFVVertex::doinit();
  tUEDBasePtr UEDBase = 
    dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!UEDBase)
    throw InitException() << "UEDF1F1P0Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theCharges[0] = UEDBase->ee();
  theCharges[1] = UEDBase->ed();
  theCharges[2] = UEDBase->eu();
}

NoPIOClassDescription<UEDF1F1P0Vertex> UEDF1F1P0Vertex::initUEDF1F1P0Vertex;
// Definition of the static class description member.

void UEDF1F1P0Vertex::Init() {

  static ClassDocumentation<UEDF1F1P0Vertex> documentation
    ("This class couples a pair of level-1 KK fermions to an SM "
     "photon.");

}


void UEDF1F1P0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long iferm;
  if(part1->id() == ParticleID::gamma) 
    iferm  = abs(part2->id());
  else if(part2->id() == ParticleID::gamma)
    iferm = abs(part1->id());
  else if(part3->id() == ParticleID::gamma)
    iferm = abs(part1->id());
  else {
    throw HelicityLogicalError() << "UEDF1F1P0Vertex::setCoupling - There is no "
				 << "photon in this vertex!"
				 << Exception::warning;
    norm(0.0);
    return;
  }
  if((iferm >= 5100001 && iferm <= 5100006) || 
     (iferm >= 5100011 && iferm <= 5100016) ||
     (iferm >= 6100001 && iferm <= 6100006) || 
     (iferm >= 6100011 && iferm <= 6100016)) {
    if(q2 != theq2Last || theCoupLast == 0. ) {
      theq2Last = q2;
      theCoupLast = electroMagneticCoupling(q2);
    }
    norm(theCoupLast);
    if(iferm != thefermLast) {
      thefermLast = iferm;
      int smtype = (iferm > 6000000) ? iferm - 6100000 : iferm - 5100000;
      if(smtype >= 11) 
	theLRLast = theCharges[0];
      else
	theLRLast = (smtype % 2 == 0) ? theCharges[2] : theCharges[1];
    }
    left(theLRLast);
    right(theLRLast);
  }
  else {
    throw HelicityLogicalError() << "UEDF1F1P0Vertex::setCoupling - There is an "
				 << "unknown particle in this vertex " << iferm
				 << Exception::warning;
    norm(0.0);
  }
}
