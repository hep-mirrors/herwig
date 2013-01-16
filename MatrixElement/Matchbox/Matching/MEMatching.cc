// -*- C++ -*-
//
// MEMatching.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEMatching class.
//

#include "MEMatching.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/PDT/EnumParticles.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

using namespace Herwig;

MEMatching::MEMatching()
  : theBornScreening(true),
    theScreeningPower(2.0) {}

MEMatching::~MEMatching() {}

IBPtr MEMatching::clone() const {
  return new_ptr(*this);
}

IBPtr MEMatching::fullclone() const {
  return new_ptr(*this);
}

double MEMatching::channelWeight(int emitter, int emission, int spectator) const {
  // do the most simple thing for the time being; needs fixing later
  if ( realCXComb()->mePartonData()[emission]->id() == ParticleID::g ) {
    Energy2 pipk = 
      realCXComb()->meMomenta()[emitter] * realCXComb()->meMomenta()[spectator];
    Energy2 pipj = 
      realCXComb()->meMomenta()[emitter] * realCXComb()->meMomenta()[emission];
    Energy2 pjpk = 
      realCXComb()->meMomenta()[emission] * realCXComb()->meMomenta()[spectator];
    return GeV2 * pipk / ( pipj * ( pipj + pjpk ) );
  }
  return
    GeV2 / (realCXComb()->meMomenta()[emitter] * realCXComb()->meMomenta()[emission]);
}

double MEMatching::channelWeight() const {
  double currentChannel = channelWeight(dipole()->realEmitter(),
					dipole()->realEmission(),
					dipole()->realSpectator());
  if ( currentChannel == 0. )
    return 0.;
  double sum = 0.;
  for ( vector<Ptr<SubtractionDipole>::ptr>::const_iterator dip =
	  dipole()->partnerDipoles().begin();
	dip != dipole()->partnerDipoles().end(); ++dip )
    sum += channelWeight((**dip).realEmitter(),
			 (**dip).realEmission(),
			 (**dip).realSpectator());
  assert(sum > 0.0);
  return currentChannel / sum;
}

double MEMatching::screeningME2() const {
  return
    pow(sqr(dipole()->lastPt())/bornXComb()->lastSHat(),screeningPower()) *
    dipole()->underlyingBornME()->me2Norm();
}

CrossSection MEMatching::dSigHatDR() const {
  double pdfFactor = 1.;
  double bornPDF = bornPDFWeight(showerScalesInSubtraction());
  double bornPDFHard = bornPDF;
  if ( showerScalesInSubtraction() )
    bornPDFHard = bornPDFWeight(false);
  if ( bornScreening() ) {
    double bornME2 = dipole()->underlyingBornME()->me2();
    double screenME2 = screeningME2();
    pdfFactor = bornME2 * bornPDFHard / ( bornME2 * bornPDF + screenME2 );
  } else {
    pdfFactor = bornPDFHard / bornPDF;
  }
  assert(realXComb()->lastME2() > 0.0);
  return
    sqr(hbarc) * 
    realXComb()->jacobian() * 
    realPDFWeight(showerScalesInSubtraction()) *
    couplingWeight(showerScalesInSubtraction()) *
    pdfFactor *
    channelWeight() * realXComb()->lastME2() /
    (2. * realXComb()->lastSHat());
}

double MEMatching::me2() const {
  double bornPDF = bornPDFWeight(showerScalesInSplitting());
  double realPDF = realPDFWeight(showerScalesInSplitting());
  assert(bornXComb()->lastME2() > 0.0);
  double den = 
    bornXComb()->lastME2() * bornPDF;
  if ( bornScreening() )
    den += screeningME2();
  double num =
    dipole()->realEmissionME()->me2() * realPDF;
  num *= pow(bornXComb()->lastSHat()/realXComb()->lastSHat(),2.*(realCXComb()->mePartonData().size())-8.);
  return 
    (num/den) *
    (bornXComb()->lastSHat()/realXComb()->lastSHat()) * 
    couplingWeight(showerScalesInSplitting());
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MEMatching::persistentOutput(PersistentOStream & os) const {
  os << theBornScreening << theScreeningPower;
}

void MEMatching::persistentInput(PersistentIStream & is, int) {
  is >> theBornScreening >> theScreeningPower;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MEMatching,Herwig::ShowerApproximation>
  describeHerwigMEMatching("Herwig::MEMatching", "HwMatchbox.so");

void MEMatching::Init() {

  static ClassDocumentation<MEMatching> documentation
    ("MEMatching implements NLO matching with matrix element correction (aka Powheg).");

  static Switch<MEMatching,bool> interfaceBornScreening
    ("BornScreening",
     "Switch on or off Born screening",
     &MEMatching::theBornScreening, true, false, false);
  static SwitchOption interfaceBornScreeningOn
    (interfaceBornScreening,
     "On",
     "Perform Born screening",
     true);
  static SwitchOption interfaceBornScreeningOff
    (interfaceBornScreening,
     "Off",
     "Do not perform Born screening",
     false);

  static Parameter<MEMatching,double> interfaceScreeningPower
    ("ScreeningPower",
     "Set the power of pt used in the screening term",
     &MEMatching::theScreeningPower, 2.0, 1.0, 0,
     false, false, Interface::lowerlim);

}

