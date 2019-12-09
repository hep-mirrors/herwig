// -*- C++ -*-
//
// UEDF1F0H1Vertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDF1F0H1Vertex class.
//

#include "UEDF1F0H1Vertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDF1F0H1Vertex::UEDF1F0H1Vertex() : theRadius(ZERO), theMw(ZERO), 
				     theSinThetaW(0.), theq2Last(ZERO),
				     theCoupLast(0.), theLeftLast(0.),
				     theRightLast(0.), theAntiLast(0),
				     theFermLast(0), theHLast(0) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

void UEDF1F0H1Vertex::doinit() {
  long heavy[3] = {5, 6, 15};
  //h0
  for( int i = 0; i < 3; ++i ) {
    addToList(-5100000 - i, 5100000 + i, 25);
    addToList(-6100000 - i, 6100000 + i, 25);
    addToList(-5100000 - i, 6100000 + i, 25);
    addToList(-6100000 - i, 5100000 + i, 25);
  }
  // Neutral KK-Higgs
  long higgs[2] = {5100025, 5100036};
  for( unsigned int h = 0; h < 2; ++h ) {
    for( int i = 0; i < 3; ++i ) {
      addToList(-heavy[i], 5100000 + heavy[i], higgs[h]);
      addToList(-5100000 - heavy[i], heavy[i], higgs[h]);
      addToList(-heavy[i], 6100000 + heavy[i], higgs[h]);
      addToList(-6100000 - heavy[i], heavy[i], higgs[h]);
    }
  }

  //KK-charged higgs
  //outgoing H+
  addToList(-5100006, 5, 5100037);
  addToList(-6100006, 5, 5100037);

  addToList(-6, 5100005, 5100037);
  addToList(-6, 6100005, 5100037);

  addToList(-5100016, 15, 5100037);
  addToList(-6100016, 15, 5100037);

  addToList(-16, 5100015, 5100037);
  addToList(-16, 6100015, 5100037);

  //outgoing H-
  addToList(-5100005, 6,-5100037);
  addToList(-6100005, 6,-5100037);

  addToList(-5, 5100006,-5100037);
  addToList(-5, 6100006,-5100037);

  addToList(-5100015, 16,-5100037);
  addToList(-6100015, 16,-5100037);

  addToList(-15, 5100016,-5100037);
  addToList(-15, 6100016,-5100037);
  FFSVertex::doinit();
  tUEDBasePtr UEDBase = 
    dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!UEDBase)
    throw InitException() << "UEDF1F0H1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theRadius = UEDBase->compactRadius();
  theSinThetaW = sqrt(sin2ThetaW());
  theCosThetaW = sqrt(1. - sin2ThetaW());
  theMw = getParticleData(24)->mass();
  theMz = getParticleData(23)->mass();
}


void UEDF1F0H1Vertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theRadius,1/GeV) << ounit(theMw,GeV) << theSinThetaW 
     << ounit(theMz, GeV) << theCosThetaW;
}

void UEDF1F0H1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theRadius,1/GeV) >> iunit(theMw,GeV) >> theSinThetaW
     >> iunit(theMz, GeV) >> theCosThetaW;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<UEDF1F0H1Vertex,FFSVertex>
describeHerwigUEDF1F0H1Vertex("Herwig::UEDF1F0H1Vertex", "HwUED.so");

void UEDF1F0H1Vertex::Init() {

  static ClassDocumentation<UEDF1F0H1Vertex> documentation
    ("The coupling involving a KK-Higgs and a pair of fermions.");

}

void UEDF1F0H1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long anti(abs(part1->id())), ferm(abs(part2->id())), higgs(part3->id());
  if( ferm > 17 ) swap( ferm, anti);

  if( anti != theAntiLast || ferm != theFermLast || higgs != theHLast ) { 
    theAntiLast = anti;
    theFermLast = ferm;
    theHLast = higgs;

    tcPDPtr pd;
    if( higgs != 25 ) {
      pd = getParticleData(ferm);
    }
    else {
      long smid = ( ferm/1000000 == 5 ) ? ferm - 5100000 : ferm - 6100000;
      pd = getParticleData(smid);
    }
    Energy mf = pd->mass();
    double alpha = mf*theRadius/2.;
    double salpha = sin(alpha);
    double calpha = cos(alpha);
    double fact(0.);
    if( abs(higgs) == 5100037 ) {
      fact = theRadius/2./sqrt(1. + sqr(theMw*theRadius)) * UnitRemoval::E;
      theRightLast = theRadius*theMw;
      if(anti/1000000 == 5) {
	Energy mfk  = getParticleData(anti - 5100000)->mass();
	theLeftLast = (theMw*calpha - (mfk*salpha/theRadius/theMw)) 
	  * UnitRemoval::InvE;
      
	theRightLast *= calpha*mfk* UnitRemoval::InvE;
      }
      else {
	Energy mfk = getParticleData(anti - 6100000)->mass();
	theLeftLast = (theMw*salpha +(mfk*calpha/theRadius/theMw))
	  * UnitRemoval::InvE;
	theRightLast *= -salpha*mfk*UnitRemoval::InvE;
      }
      theLeftLast *= fact;
      theRightLast *= fact;
      if( higgs < 0 ) swap( theLeftLast, theRightLast );
    }
    else if( higgs == 5100025 ) {
      fact = mf/theMw/2.;
      if( anti/1000000 == 5 )
	theLeftLast = salpha + calpha;
      else 
	theLeftLast = salpha - calpha;
      theRightLast = theLeftLast;

      theLeftLast *= fact;
      theRightLast *= fact;
    }
    else if( higgs == 5100036 ) {
      fact = theRadius/theCosThetaW/sqrt(1.+sqr(theMz*theRadius))*UnitRemoval::E;
    
      double i3f = ( ferm % 2 == 0 ) ? 0.5 : -0.5;
      double qf = pd->charge()/eplus;
      if( anti/1000000 == 5 ) {
	theLeftLast = (theMz*calpha*(i3f - qf*sqr(theSinThetaW))
		       - mf*salpha/theRadius/theMw) * UnitRemoval::InvE;
	theRightLast = (-theMz*salpha*qf*sqr(theSinThetaW) 
			+ mf*calpha/theRadius/theMw) * UnitRemoval::InvE; 
      }
      else {
	theLeftLast = (theMz*salpha*(i3f - qf*sqr(theSinThetaW))
		       - mf*calpha/theRadius/theMw) * UnitRemoval::InvE;
	theRightLast = (-theMz*calpha*qf*sqr(theSinThetaW) 
			+ mf*salpha/theRadius/theMw)*UnitRemoval::InvE; 
      }
      theLeftLast *= fact;
      theRightLast *= fact;
    }
    else {
      theLeftLast = mf*calpha*salpha/2./theMw;
      if( ferm/1000000 == 5 ) theLeftLast *= -1.;  
      theRightLast = theLeftLast;
    }
  }


  if(q2 != theq2Last || theCoupLast == 0.) {
    theq2Last = q2;
    theCoupLast = weakCoupling(q2);
  }

  norm(theCoupLast);
  left(theLeftLast);
  right(theRightLast);
}
