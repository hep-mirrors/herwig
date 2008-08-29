// -*- C++ -*-
//
// UEDF1F1Z0Vertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F1Z0Vertex class.
//

#include "UEDF1F1Z0Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F1Z0Vertex::UEDF1F1Z0Vertex() : theSin2ThW(0.0), theCosThW(0.0), theRadius(),
				     theID1Last(0), theID2Last(0) ,
				     theq2Last(0.*GeV2), theCoupLast(0.), 
				     theLeftLast(0.), theRightLast(0.) {
  vector<long> anti, ferm, boson(25, 23);
  //QQ, uu, dd
  for(long i = 5100001; i < 6100007; ++i) {
    if(i == 5100007) i += 999994;
    anti.push_back(-i);
    ferm.push_back(i);
  }
  //top/bottom quark l/r mixing
  anti.push_back(-5100006); ferm.push_back(6100006); 
  anti.push_back(-6100006); ferm.push_back(5100006); 
  anti.push_back(-5100005); ferm.push_back(6100005); 
  anti.push_back(-6100005); ferm.push_back(5100005); 
  //leptons
  for(long i = 5100011; i < 5100017; ++i) {
    anti.push_back(-i);
    ferm.push_back(i);
  }
  for(long i = 6100011; i < 6100017; i +=2) {
    anti.push_back(-i);
    ferm.push_back(i);
  }
  setList(anti, ferm, boson);
}

void UEDF1F1Z0Vertex::doinit() throw(InitException) {
  FFVVertex::doinit();
  UEDBasePtr model = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "UEDF1F1Z0Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  
  theSin2ThW = model->sin2ThetaW();
  theCosThW = sqrt(1. - theSin2ThW); 
  theRadius = model->compactRadius();
  orderInGs(0);
  orderInGem(1);
}

void UEDF1F1Z0Vertex::persistentOutput(PersistentOStream & os) const {
  os << theSin2ThW << theCosThW << ounit(theRadius,1/GeV);
}

void UEDF1F1Z0Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theSin2ThW >> theCosThW >> iunit(theRadius,1/GeV);
}

ClassDescription<UEDF1F1Z0Vertex> UEDF1F1Z0Vertex::initUEDF1F1Z0Vertex;
// Definition of the static class description member.

void UEDF1F1Z0Vertex::Init() {

  static ClassDocumentation<UEDF1F1Z0Vertex> documentation
    ("This is the implementation of the level-1 fermion pair Z_0 boson "
     "coupling.");

}

void UEDF1F1Z0Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  if( part3->id() != 23 ) {
    setNorm(0.0);
    setLeft(0.0);
    setRight(0.0);  
    throw HelicityLogicalError()
      << "UEDF1F1Z0Vertex::setCoupling - The vector boson in this vertex "
      << "is not a Z^0 boson. ID: " << part3->id() << "\n"
      << Exception::warning;
    return;
  }
  long ianti(abs(part1->id())), iferm(abs(part2->id()));
  bool ferma = (iferm >= 5100001 && iferm <= 5100006) ||
    (iferm >= 6100001 && iferm <= 6100006) || 
    (iferm >= 5100011 && iferm <= 5100016) ||
    (iferm >= 6100011 && iferm <= 6100016); 
  bool fermb = (ianti >= 5100001 && ianti <= 5100006) ||
    (ianti >= 6100001 && ianti <= 6100006) || 
    (ianti >= 5100011 && ianti <= 5100016) ||
    (ianti >= 6100011 && ianti <= 6100016);
  if( ferma && fermb  ) {
    if(q2 != theq2Last) {
	theq2Last = q2;
	theCoupLast = 0.5*weakCoupling(q2)/theCosThW;
    }
    if( ianti != theID1Last || iferm != theID2Last) {
      theID1Last = ianti;
      theID2Last = iferm;
      int stateA = ianti/1000000;
      int stateB = iferm/1000000;
      long smID = (stateA == 6) ? ianti - 6100000 : ianti - 5100000;
      // L/R mixing
      double alpha = atan(getParticleData(smID)->mass()*theRadius)/2.;
      double sin2al = sqr(sin(alpha));
      double cos2al = 1. - sin2al;
      
      if(stateA == 5 && stateB == 5) {
	if(smID >= 11 && smID <= 16)
	  theLeftLast = -cos2al + 2.*theSin2ThW;
	else if(smID <= 6 && smID % 2 == 0)
	  theLeftLast = cos2al - 4.*theSin2ThW/3.;
	else
	  theLeftLast = -cos2al + 2.*theSin2ThW/3.;

	theRightLast = theLeftLast;
      }
      else if(stateA == 6 && stateB == 6) {
	if(smID >= 11 && smID <= 16)
	  theLeftLast = -sin2al + 2.*theSin2ThW;
	else if(smID <=6 && smID % 2 == 0)
	  theLeftLast = sin2al - 4.*theSin2ThW/3.;
	else
	  theLeftLast = -sin2al + 2.*theSin2ThW/3.;

	theRightLast = theLeftLast;
      }
      else {
	theLeftLast = sqrt(sin2al*cos2al);
	if(smID % 2 == 0) theLeftLast *= -1.;
	theRightLast = -theLeftLast;
      }
    }
    setNorm(theCoupLast);
    setLeft(theLeftLast);
    setRight(theRightLast);
  }
  else {
    throw HelicityLogicalError() << "UEDF1F1Z0Vertex::setCoupling - "
				 << "There is an unknown particle(s) in the "
				 << "UED F^(1) F^(1) Z^(0) vertex. ID: " 
				 << ianti << " " << iferm 
				 << Exception::warning;      
    setNorm(0.0);
    setLeft(0.0);
    setRight(0.0);  
  }
}
