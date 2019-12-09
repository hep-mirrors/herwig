// -*- C++ -*-
//
// DipoleMatching.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleMatching class.
//

#include "DipoleMatching.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

using namespace Herwig;

DipoleMatching::DipoleMatching() {}

DipoleMatching::~DipoleMatching() {}

IBPtr DipoleMatching::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleMatching::fullclone() const {
  return new_ptr(*this);
}

CrossSection DipoleMatching::dSigHatDR() const {

  double xme2 = 0.;

  pair<int,int> ij(dipole()->bornEmitter(),
		   dipole()->bornSpectator());
  double ccme2 = 
    dipole()->underlyingBornME()->largeNColourCorrelatedME2(ij,theLargeNBasis);

   if(ccme2==0.)return 0.*nanobarn;
      
   double lnme2=dipole()->underlyingBornME()->largeNME2(theLargeNBasis);
   if(lnme2==0){
     generator()->log() <<"\nDipoleMatching: ";
     generator()->log() <<"\n  LargeNME2 is ZERO, while largeNColourCorrelatedME2 is not ZERO." ;
     generator()->log() <<"\n  This is too seriuos.\n" ;
     generator()->log() << Exception::runerror;
   }    
    
    
  ccme2 *=
    dipole()->underlyingBornME()->me2() /lnme2;

  xme2 = dipole()->me2Avg(ccme2);

  xme2 /= dipole()->underlyingBornME()->lastXComb().lastAlphaS();
  double bornPDF = bornPDFWeight(dipole()->underlyingBornME()->lastScale());
  if ( bornPDF == 0.0 )
    return ZERO;
  xme2 *= bornPDF;

  if ( profileScales() )
    xme2 *= profileScales()->hardScaleProfile(dipole()->showerHardScale(),dipole()->lastPt());

  CrossSection res = 
    sqr(hbarc) * 
    realXComb()->jacobian() * 
    subtractionScaleWeight() *
    xme2 /
    (2. * realXComb()->lastSHat());

  return res;

}

double DipoleMatching::me2() const {
  throw Exception() << "DipoleMatching::me2(): Not intented to use. Disable the ShowerApproximationGenerator."
		    << Exception::runerror;
  return 0.;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleMatching::persistentOutput(PersistentOStream & os) const {
  os << theShowerHandler;
}

void DipoleMatching::persistentInput(PersistentIStream & is, int) {
  is >> theShowerHandler;
}

void DipoleMatching::doinit() {
  if ( theShowerHandler ) {
    theShowerHandler->init();
    hardScaleFactor(theShowerHandler->hardScaleFactor());
    factorizationScaleFactor(theShowerHandler->factorizationScaleFactor());
    renormalizationScaleFactor(theShowerHandler->renormalizationScaleFactor());
    profileScales(theShowerHandler->profileScales());
    restrictPhasespace(theShowerHandler->restrictPhasespace());
    hardScaleIsMuF(theShowerHandler->hardScaleIsMuF());
    if ( theShowerHandler->showerPhaseSpaceOption() == 0 ) {
      useOpenZ(false);
    } else if ( theShowerHandler->showerPhaseSpaceOption() == 1 ) {
      useOpenZ(true);
    } else {
      throw InitException() << "DipoleMatching::doinit(): Choice of shower phase space cannot be handled by the matching";
    }
  }
  // need to fo this after for consistency checks
  ShowerApproximation::doinit();
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<DipoleMatching,Herwig::ShowerApproximation>
  describeHerwigDipoleMatching("Herwig::DipoleMatching", "HwDipoleMatching.so HwShower.so");

void DipoleMatching::Init() {

  static ClassDocumentation<DipoleMatching> documentation
    ("DipoleMatching implements NLO matching with the dipole shower.");

  static Reference<DipoleMatching,ShowerHandler> interfaceShowerHandler
    ("ShowerHandler",
     "The dipole shower handler object to use.",
     &DipoleMatching::theShowerHandler, false, false, true, true, false);
  interfaceShowerHandler.rank(-1);

}

