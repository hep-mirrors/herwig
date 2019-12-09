// -*- C++ -*-
//
// UEDW0W1W1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDW0W1W1Vertex class.
//

#include "UEDW0W1W1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDW0W1W1Vertex::UEDW0W1W1Vertex() : theSinW(0.), theCosW(0.),
				     theSinThetaOne(0.), theCosThetaOne(0.),
				     theq2last(), theElast(0.), theCouplast(0.),
				     theSMlast(0), theKKlast(0) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void UEDW0W1W1Vertex::doinit() {
  addToList( 22, -5100024, 5100024);
  addToList( 23, -5100024, 5100024);

  addToList( 24, -5100024, 5100022);
  addToList( 24, -5100024, 5100023);

  addToList(-24,  5100024, 5100022);
  addToList(-24,  5100024, 5100023);
  VVVVertex::doinit();
  tUEDBasePtr model = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "UEDW0W1W1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theSinW = sqrt(sin2ThetaW());
  theCosW = sqrt( 1. - sqr(theSinW) );
  theSinThetaOne = model->sinThetaOne();
  theCosThetaOne = sqrt( 1. - sqr(theSinThetaOne));
}

void UEDW0W1W1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theSinW << theCosW << theSinThetaOne << theCosThetaOne;
}

void UEDW0W1W1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theSinW >> theCosW >> theSinThetaOne >> theCosThetaOne;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDW0W1W1Vertex,VVVVertex>
describeHerwigUEDW0W1W1Vertex("Herwig::UEDW0W1W1Vertex", "HwUED.so");

void UEDW0W1W1Vertex::Init() {

  static ClassDocumentation<UEDW0W1W1Vertex> documentation
    ("The coupling of an SM W boson to a level 1 KK W and KK Z and KK photon");

}

/// \todo look again
void UEDW0W1W1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long id1(abs(part1->id())), id2(abs(part2->id())), id3(abs(part3->id())), 
    smID(0), kkparticle(0);
  double perm(1.);
  if( id1 == 22 || id1 == 23) {
    smID = id1;
    kkparticle = id2;
    if(part2->id()>0) perm=-1.;
  }
  else if(id2 == 22 || id2 == 23) {
    smID = id2;
    kkparticle = id1;
    if(part1->id()<0) perm=-1.;
  }
  else if(id3 == 22 || id3 == 23) {
    smID = id3;
    kkparticle = id1;
    if(part1->id()>0) perm=-1.;
  }
  else if(id1 == 24 ) {
    if( part1->id() == 24 ) perm = -1.;
    smID = id1;
    kkparticle = (id2 == 5100024) ? id3 : id2;
    if(id3 == 5100024) perm *=-1.;
  }
  else if( id2 == 24 ) {
    if( part2->id() == 24 ) perm = -1.;
    smID = id2;
    kkparticle = (id1 == 5100024) ? id3 : id1;
    if(id1 == 5100024) perm *=-1.;
  }
  else if( id3 == 24 ) {
    if( part3->id() == 24 ) perm = -1.;
    smID = id3;
    kkparticle = (id1 == 5100024) ? id2 : id1;
    if(id2 == 5100024) perm *=-1.;
  }
  else {
    throw HelicityLogicalError()
      << "UEDW0W1W1Vertex::setCoupling() - There is no SM gauge boson in "
      << "this vertex. " << id1 << " " << id2 << " " << id3 
      << Exception::warning; 
    norm(0.);
    return;
  }
  if( q2 != theq2last || theElast == 0.) {
    theq2last = q2;
    theElast = electroMagneticCoupling(q2);
  }
  
  if( smID != theSMlast || kkparticle != theKKlast ) { 
    theSMlast = smID;
    theKKlast = kkparticle;
    if( smID == 22 )
      theCouplast = 1.;
    else if(smID == 23) 
      theCouplast = theCosW/theSinW;
    else {
      if( kkparticle == 5100023 )
	theCouplast = theCosThetaOne/theSinW;
      else
	theCouplast = theSinThetaOne/theSinW;
    }
  }
  norm(perm*theElast*theCouplast);
}
