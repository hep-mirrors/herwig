// -*- C++ -*-
//
// UEDF1F1P0Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F1P0Vertex class.
//

#include "UEDF1F1P0Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F1P0Vertex::UEDF1F1P0Vertex() : coupLast_(0.0), q2Last_(ZERO),
				     fermLast_(0), LRLast_(0.0), 
				     charges_(3) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

void UEDF1F1P0Vertex::persistentOutput(PersistentOStream & os) const {
  os << charges_;
}

void UEDF1F1P0Vertex::persistentInput(PersistentIStream & is, int) {
  is >> charges_;
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
  charges_[0] = UEDBase->ee();
  charges_[1] = UEDBase->ed();
  charges_[2] = UEDBase->eu();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDF1F1P0Vertex,FFVVertex>
describeHerwigUEDF1F1P0Vertex("Herwig::UEDF1F1P0Vertex", "HwUED.so");

void UEDF1F1P0Vertex::Init() {

  static ClassDocumentation<UEDF1F1P0Vertex> documentation
    ("This class couples a pair of level-1 KK fermions to an SM "
     "photon.");

}


void UEDF1F1P0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr ,
#ifndef NDEBUG
				  tcPDPtr part3) {
#else
				  tcPDPtr ) {
#endif
  long iferm = abs(part1->id());
  assert(part3->id()==ParticleID::gamma);
  assert((iferm >= 5100001 && iferm <= 5100006) || 
	 (iferm >= 5100011 && iferm <= 5100016) ||
	 (iferm >= 6100001 && iferm <= 6100006) || 
	 (iferm >= 6100011 && iferm <= 6100016));
  if(q2 != q2Last_ || coupLast_ == 0. ) {
    q2Last_ = q2;
    coupLast_ = -electroMagneticCoupling(q2);
  }
  norm(coupLast_);
  if(iferm != fermLast_) {
    fermLast_ = iferm;
    int smtype = (iferm > 6000000) ? iferm - 6100000 : iferm - 5100000;
    if(smtype >= 11) 
      LRLast_ = charges_[0];
    else
      LRLast_ = (smtype % 2 == 0) ? charges_[2] : charges_[1];
  }
  left(LRLast_);
  right(LRLast_);
}
