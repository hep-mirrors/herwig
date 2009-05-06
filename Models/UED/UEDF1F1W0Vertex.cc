// -*- C++ -*-
//
// UEDF1F1W0Vertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F1W0Vertex class.
//

#include "UEDF1F1W0Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F1W0Vertex::UEDF1F1W0Vertex(): theRadius(ZERO), theQ2Last(ZERO), 
				    theCoupLast(0.), 
				    thefermALast(0), thefermBLast(0) {
  vector<long> anti, ferm, wboson;
  //outgoing W+
  for( long i = 2; i < 17; i += 2 ) {
    if( i == 7 ) i += 5;
    anti.push_back(-5000000 - i);
    ferm.push_back(5000000 + i - 1);
    wboson.push_back(24);
    if( i < 7 ) {
      anti.push_back(-6000000 - i);
      ferm.push_back(6000000 + i - 1);
      wboson.push_back(24);
    }
  }
  anti.push_back(-6000006);
  ferm.push_back(5000005);
  wboson.push_back(24);
  anti.push_back(-5000005);
  ferm.push_back(6000006);
  wboson.push_back(24);
  //outgoing W-
  for( long i = 1; i < 16; i += 2 ) {
    if( i == 6 ) i += 5;
    anti.push_back(-5000000 - i);
    ferm.push_back(5000001 + i);
    wboson.push_back(-24);
    if( i < 6 ) {
      anti.push_back(-6000000 - i);
      ferm.push_back(6000001 + i);
      wboson.push_back(-24);
    }
  }
  anti.push_back(-6000005);
  ferm.push_back(5000006);
  wboson.push_back(-24);
  anti.push_back(-5000005);
  ferm.push_back(6000006);
  wboson.push_back(-24);
  
  setList(anti, ferm, wboson);
}

void UEDF1F1W0Vertex::doinit() {
  FFVVertex::doinit();
  tUEDBasePtr model = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "UEDF1F1W0Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theRadius = model->compactRadius();
  orderInGs(0);
  orderInGem(0);
}

void UEDF1F1W0Vertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theRadius,1/GeV);
}

void UEDF1F1W0Vertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theRadius,1/GeV);
}

ClassDescription<UEDF1F1W0Vertex> UEDF1F1W0Vertex::initUEDF1F1W0Vertex;
// Definition of the static class description member.

void UEDF1F1W0Vertex::Init() {

  static ClassDocumentation<UEDF1F1W0Vertex> documentation
    ("This class implements the coupling of a pair of level-1 KK fermions"
     "to an SM W boson");

}

void UEDF1F1W0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
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
  if(ferma && fermb) {
    if(q2 != theQ2Last || theCoupLast == 0. ) {
      theQ2Last = q2;
      theCoupLast = sqrt(0.5)*weakCoupling(q2);
    }
    if(iferm != thefermALast || ianti != thefermBLast) {
      thefermALast = iferm;
      thefermBLast = ianti;
      int stateA(ianti/1000000), stateB(iferm/1000000);
      long sma = (stateA == 6) ? iferm - 6100000 : iferm - 5100000;
      long smb = (stateB == 6) ? ianti - 6100000 : ianti - 5100000;
      double afu(0.), afd(0.);
      if( sma % 2 == 0 ) {
	afu = atan(getParticleData(sma)->mass()*theRadius)/2.;
	afd = atan(getParticleData(smb)->mass()*theRadius)/2.;
      }
      else {
	afd = atan(getParticleData(sma)->mass()*theRadius)/2.;
	afu = atan(getParticleData(smb)->mass()*theRadius)/2.;
	
      }
      if( stateA == stateB ) {
	if( stateA == 5 ) {
	  setLeft(cos(afu)*cos(afd));
	  setRight(cos(afu)*cos(afd));
	}
	else {
	  setLeft(sin(afu)*sin(afd));
	  setRight(sin(afu)*sin(afd));
	}
      }
      else {
	if( sma % 2 == 0 ) {
	  if( stateA == 5 ) {
	    setLeft(cos(afu)*sin(afd));
	    setRight(-cos(afu)*sin(afd));
	  }
	  else {
	    setLeft(sin(afu)*cos(afd));
	    setRight(-sin(afu)*cos(afd));
	  }
	}
	else {
	  if( stateA == 5 ) {
	    setLeft(sin(afu)*cos(afd));
	    setRight(-sin(afu)*cos(afd));
	  }
	  else {
	    setLeft(cos(afu)*sin(afd));
	    setRight(-cos(afu)*sin(afd));
	  }
	}
      }
    }
    setNorm(theCoupLast);
  }
  else {
    throw HelicityLogicalError() << "UEDF1F1W0Vertex::setCoupling - "
				 << "There is an unknown particle(s) in the "
				 << "UED F^(1) F^(1) W^(0) vertex. ID: " 
				 << ianti << " " << iferm 
				 << Exception::warning;      
    setNorm(0.0);
    setLeft(0.0);
    setRight(0.0);  
  }
}
