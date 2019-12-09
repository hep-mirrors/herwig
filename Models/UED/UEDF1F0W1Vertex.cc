// -*- C++ -*-
//
// UEDF1F0W1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0W1Vertex class.
//

#include "UEDF1F0W1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F0W1Vertex::UEDF1F0W1Vertex() : theSinW(0.), theCosW(0.), theSinOne(0.),
				     theCosOne(0.), theSinWmO(0.), 
				     theCosWmO(0.), 
				     theCKM(0, vector<Complex>(0, 0.)),
				     theq2last(),
				     theCouplast(0.), theLlast(0.),
				     theRlast(0.), theGBlast(0), 
				     theKKlast(0), theSMlast(0) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

void UEDF1F0W1Vertex::doinit() {
  //outgoing W+
  for(long i = 2; i < 7; i += 2) {
    for(long j = 1; j < 6; j += 2) {
      addToList( -i, 5100000 + j, 5100024 );
      addToList( -(5100000 + i), j, 5100024 );
    }
  }
  for(long i = 11; i < 17; i += 2) {
    addToList( -i-1, 5100000 + i, 5100024 );
    addToList( -(5100001 + i), i, 5100024 );
  }
  //outgoing W-
  for(long i = 1; i < 6; i += 2) {
    for(long j = 2 ; j < 7; j += 2) {
      addToList( -i, 5100000 + j, -5100024 );
      addToList( -(5100000 + i), j, -5100024 );
    }
  }
  for(long i = 11; i < 17; i += 2) {
    addToList( -i, 5100001 + i, -5100024 );
    addToList(-(5100000 + i), i + 1, -5100024);
  }
  long boson[2] = {5100022,5100023}; 
  for(long b = 0; b < 2; ++b) { 
    //QQ
    for(int i = 1; i < 7; ++i) {
      addToList( -i, i + 5100000, boson[b]);
      addToList(-(i + 5100000), i, boson[b]);

      addToList(-i, i + 6100000, boson[b]);
      addToList(-(i + 6100000), i, boson[b]);
    }
    //LL
    for(int i = 11; i < 17; ++i) {
      addToList( -i, i + 5100000, boson[b]);
      addToList(-(i + 5100000), i, boson[b]);
      if( i % 2 != 0 ) {
	addToList(-i, i + 6100000, boson[b]);
	addToList(-(i + 6100000), i, boson[b]);
      }
    }
  }
  FFVVertex::doinit();
  tUEDBasePtr UEDBase = 
    dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!UEDBase)
    throw InitException() << "UEDF1F0W1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theSinW = sqrt(sin2ThetaW());
  theCosW = sqrt( 1. - sqr(theSinW));
  theSinOne = UEDBase->sinThetaOne();
  theCosOne = sqrt(1. - sqr(theSinOne)); 
  theSinWmO = theSinW*theCosOne - theSinOne*theCosW;
  theCosWmO = theCosW*theCosOne + theSinW*theSinOne;
  theCKM = dynamic_ptr_cast<Ptr<StandardCKM>::transient_pointer>
    (UEDBase->CKM())->getUnsquaredMatrix(UEDBase->families());
}


void UEDF1F0W1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theSinW << theCosW << theSinOne << theCosOne
     << theSinWmO << theCosWmO << theCKM;
}

void UEDF1F0W1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theSinW >> theCosW >> theSinOne >> theCosOne
     >> theSinWmO >> theCosWmO >> theCKM;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDF1F0W1Vertex,FFVVertex>
describeHerwigUEDF1F0W1Vertex("Herwig::UEDF1F0W1Vertex", "HwUED.so");

void UEDF1F0W1Vertex::Init() {

  static ClassDocumentation<UEDF1F0W1Vertex> documentation
    ("This is the coupling of a KK1 W boson to a KK1 fermion and "
     "a SM fermion.");

}

void UEDF1F0W1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long id1(abs(part1->id())), id2(abs(part2->id())),
    gboson(abs(part3->id())), kkparticle(0), smID(0);
  assert( gboson == 5100022  || gboson == 5100023 || gboson == 5100024 );
  if( id1 > 5000000 ) {
    kkparticle = id1;
    smID = id2;
  }
  else {
    kkparticle = id2;
    smID = id1;
  }
  if( (kkparticle >= 5100001 && kkparticle <= 5100006) ||
      (kkparticle >= 6100001 && kkparticle <= 6100006) ||
      (kkparticle >= 5100011 && kkparticle <= 5100016) ||
      (kkparticle >= 6100011 && kkparticle <= 6100016) ) {
    if(q2 != theq2last || theCouplast == 0.) {
      theq2last = q2;
      theCouplast = electroMagneticCoupling(q2);
    }
    if( gboson != theGBlast || kkparticle != theKKlast || smID != theSMlast ) {
      theGBlast = gboson;
      theKKlast = kkparticle;
      theSMlast = smID;
      if( gboson == 5100024 ) {
	Complex ckm(1.);
	if( smID >= 1 && smID <= 6 ) {
	  long smIDb(kkparticle - 5100000);
	  if( smID % 2 != 0 ) swap(smID, smIDb);
	  ckm = theCKM[smID/2 - 1][(smIDb - 1)/2];
	}
	theLlast = -ckm/sqrt(2)/theSinW;
	theRlast = 0.;
      }
      else if( gboson == 5100022 || gboson == 5100023 ) {
	double Qf = getParticleData(smID)->charge()/eplus;
	if( kkparticle/1000000 == 5 ) {
	  theRlast = 0.;
	  double I3f = (abs(smID) % 2 == 0) ? 0.5 : -0.5;
	  if( gboson == 5100023 )
	    theLlast = (Qf*theSinOne 
			- I3f*theCosWmO/theSinW)/theCosW;
	  else
	    theLlast = -(Qf*theCosOne 
			 - I3f*theSinWmO/theSinW)/theCosW;
	}
	else {
	  theLlast = 0.;
	  if( gboson == 5100023 )
	    theRlast = Qf*theSinOne/theCosW;
	  else
	    theRlast = -Qf*theCosOne/theCosW;
	}
      }
    }
    norm(theCouplast);
    left(theLlast);
    right(theRlast);
  }
  else
    throw HelicityLogicalError() << "UEDF1F0W1Vertex::setCoupling - "
				 << "There is an unknown particle in this vertex! "
				 << kkparticle
				 << Exception::warning;
}

