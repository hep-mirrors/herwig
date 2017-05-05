// -*- C++ -*-
//
// UEDF1F1W0Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F1W0Vertex class.
//

#include "UEDF1F1W0Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F1W0Vertex::UEDF1F1W0Vertex(): includeMixing_(true),
				    theRadius(ZERO), theQ2Last(ZERO), 
				    theCoupLast(0.), 
				    thefermALast(0), thefermBLast(0) {
  orderInGs(0);
  orderInGem(1);
}

void UEDF1F1W0Vertex::doinit() {
  //outgoing W+
  for( long i = 2; i < 17; i += 2 ) {
    if( i == 7 ) i += 5;
    addToList(-5100000 - i, 5100000 + i - 1, 24);
    if( i < 7 ) {
      addToList(-6100000 - i, 6100000 + i - 1, 24);
    }
  }
  if(includeMixing_) {
    addToList(-6100006, 5100005, 24);
    addToList(-5100006, 6100005, 24);
  }
  //outgoing W-
  for( long i = 1; i < 16; i += 2 ) {
    if( i == 6 ) i += 5;
    addToList(-5100000 - i, 5100001 + i, -24);
    if( i < 6 ) {
      addToList(-6100000 - i, 6100001 + i, -24);
    }
  }
  if(includeMixing_) {
    addToList(-6100005, 5100006, -24);
    addToList(-5100005, 6100006, -24);
  }
  FFVVertex::doinit();
  tUEDBasePtr model = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "UEDF1F1W0Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theRadius = model->compactRadius();
}

void UEDF1F1W0Vertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theRadius,1/GeV) << includeMixing_;
}

void UEDF1F1W0Vertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theRadius,1/GeV) >> includeMixing_;
}

ClassDescription<UEDF1F1W0Vertex> UEDF1F1W0Vertex::initUEDF1F1W0Vertex;
// Definition of the static class description member.

void UEDF1F1W0Vertex::Init() {

  static ClassDocumentation<UEDF1F1W0Vertex> documentation
    ("This class implements the coupling of a pair of level-1 KK fermions"
     "to an SM W boson");

  static Switch<UEDF1F1W0Vertex,bool> interfaceIncludeMixing
    ("IncludeMixing",
     "Include the mixing",
     &UEDF1F1W0Vertex::includeMixing_, true, false, false);
  static SwitchOption interfaceIncludeMixingYes
    (interfaceIncludeMixing,
     "Yes",
     "Include mixing",
     true);
  static SwitchOption interfaceIncludeMixingNo
    (interfaceIncludeMixing,
     "No",
     "Don't include mixing",
     false);

}

void UEDF1F1W0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
#ifndef NDEBUG
				  tcPDPtr part3) {
#else
				  tcPDPtr) {
#endif
  long ianti(abs(part1->id())), iferm(abs(part2->id()));
  assert( abs(part3->id()) == 24 );
  bool ferma = (iferm >= 5100001 && iferm <= 5100006) ||
    (iferm >= 6100001 && iferm <= 6100006) || 
    (iferm >= 5100011 && iferm <= 5100016) ||
    (iferm >= 6100011 && iferm <= 6100016); 
  bool fermb = (ianti >= 5100001 && ianti <= 5100006) ||
    (ianti >= 6100001 && ianti <= 6100006) || 
    (ianti >= 5100011 && ianti <= 5100016) ||
    (ianti >= 6100011 && ianti <= 6100016);
  if( !ferma || !fermb ) 
    throw HelicityLogicalError() << "UEDF1F1W0Vertex::setCoupling - "
				 << "There is an unknown particle(s) in the "
				 << "UED F^(1) F^(1) W^(0) vertex. ID: " 
				 << ianti << " " << iferm 
				 << Exception::runerror;
  if(q2 != theQ2Last || theCoupLast == 0. ) {
    theQ2Last = q2;
    theCoupLast = sqrt(0.5)*weakCoupling(q2);
  }
  if(iferm != thefermALast || ianti != thefermBLast) {
    thefermALast = iferm;
    thefermBLast = ianti;
    int stateA(ianti/1000000), stateB(iferm/1000000);
    long sma = (stateA == 6) ? ianti - 6100000 : ianti - 5100000;
    long smb = (stateB == 6) ? iferm - 6100000 : iferm - 5100000;
    double afu(0.), afd(0.);
    if(includeMixing_) {
      if( sma % 2 == 0 ) {
	afu = atan(getParticleData(sma)->mass()*theRadius)/2.;
	afd = atan(getParticleData(smb)->mass()*theRadius)/2.;
      }
      else {
	afd = atan(getParticleData(sma)->mass()*theRadius)/2.;
	afu = atan(getParticleData(smb)->mass()*theRadius)/2.;
      }
    }
    else {
      afd = afu = 0.;
    }
    if( stateA == stateB ) {
      if( stateA == 5 ) {
	left(cos(afu)*cos(afd));
	right(cos(afu)*cos(afd));
      }
      else {
	left(sin(afu)*sin(afd));
	right(sin(afu)*sin(afd));
      }
    }
    else {
      if( sma % 2 == 0 ) {
	if( stateA == 5 ) {
	  left(cos(afu)*sin(afd));
	  right(-cos(afu)*sin(afd));
	}
	else {
	  left(sin(afu)*cos(afd));
	  right(-sin(afu)*cos(afd));
	}
      }
      else {
	if( stateA == 5 ) {
	  left(sin(afu)*cos(afd));
	  right(-sin(afu)*cos(afd));
	}
	else {
	  left(cos(afu)*sin(afd));
	  right(-cos(afu)*sin(afd));
	}
      }
    }
  }
  norm(theCoupLast);
}
