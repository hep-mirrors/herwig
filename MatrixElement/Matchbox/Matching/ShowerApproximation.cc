// -*- C++ -*-
//
// ShowerApproximation.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerApproximation class.
//

#include "ShowerApproximation.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

using namespace Herwig;

ShowerApproximation::ShowerApproximation() 
  : HandlerBase(), theBelowCutoff(false),
    theFFPtCut(1.0*GeV), theFIPtCut(1.0*GeV), theIIPtCut(1.0*GeV),
    theShowerScalesInSubtraction(false),
    theShowerScalesInSplitting(true),
    theRestrictPhasespace(true), theHardScaleFactor(1.0),
    theExtrapolationX(0.65) {}

ShowerApproximation::~ShowerApproximation() {}

void ShowerApproximation::setDipole(Ptr<SubtractionDipole>::tcptr dip) { theDipole = dip; }

Ptr<SubtractionDipole>::tcptr ShowerApproximation::dipole() const { return theDipole; }

bool ShowerApproximation::isAboveCutoff() const {

  if ( dipole()->bornEmitter() > 1 &&
       dipole()->bornSpectator() > 1 ) {
    return dipole()->lastPt() > ffPtCut();
  } else if ( ( dipole()->bornEmitter() > 1 &&
		dipole()->bornSpectator() < 2 ) ||
	      ( dipole()->bornEmitter() < 2 &&
		dipole()->bornSpectator() > 1 ) ) {
    return dipole()->lastPt() > fiPtCut();
  } else {
    assert(dipole()->bornEmitter() < 2 &&
	   dipole()->bornSpectator() < 2);
    return dipole()->lastPt() > iiPtCut();
  }

  return true;

}

bool ShowerApproximation::isInShowerPhasespace() const {

  if ( !isAboveCutoff() )
    return false;
  if ( !restrictPhasespace() )
    return true;

  Energy maxPt = generator()->maximumCMEnergy();
  vector<Lorentz5Momentum>::const_iterator p = 
    bornCXComb()->meMomenta().begin() + 2;
  cPDVector::const_iterator pp = 
    bornCXComb()->mePartonData().begin() + 2;
  for ( ; p != bornCXComb()->meMomenta().end(); ++p, ++pp )
    if ( (**pp).coloured() )
      maxPt = min(maxPt,p->perp());
  if ( maxPt == generator()->maximumCMEnergy() )
    maxPt = (bornCXComb()->meMomenta()[0] + bornCXComb()->meMomenta()[1]).m();
  maxPt *= sqrt(hardScaleFactor());

  return dipole()->lastPt() <= maxPt;

}

double ShowerApproximation::couplingWeight(bool showerscales) const {
  if ( !showerscales )
    return 1.;
  double hardAlpha = dipole()->realEmissionME()->lastAlphaS();
  Energy2 mur = sqr(dipole()->lastPt());
  mur *= dipole()->realEmissionME()->renormalizationScaleFactor();
  double runAlpha = SM().alphaS(mur);
  return runAlpha/hardAlpha;
}

double ShowerApproximation::bornPDFWeight(bool showerscales) const {
  if ( !bornCXComb()->mePartonData()[0]->coloured() &&
       !bornCXComb()->mePartonData()[1]->coloured() )
    return 1.;
  Energy2 muf;
  if ( showerscales ) {
    muf = sqr(dipole()->lastPt());
  } else {
    muf = dipole()->underlyingBornME()->factorizationScale();
  }
  muf *= dipole()->underlyingBornME()->factorizationScaleFactor();
  double pdfweight = 1.;
  if ( bornCXComb()->mePartonData()[0]->coloured() &&
       dipole()->underlyingBornME()->havePDFWeight1() )
    pdfweight *= dipole()->underlyingBornME()->pdf1(muf,theExtrapolationX);
  if ( bornCXComb()->mePartonData()[1]->coloured() &&
       dipole()->underlyingBornME()->havePDFWeight2() )
    pdfweight *= dipole()->underlyingBornME()->pdf2(muf,theExtrapolationX);
  return pdfweight;
}

double ShowerApproximation::realPDFWeight(bool showerscales) const {
  if ( !realCXComb()->mePartonData()[0]->coloured() &&
       !realCXComb()->mePartonData()[1]->coloured() )
    return 1.;
  Energy2 muf;
  if ( showerscales ) {
    muf = sqr(dipole()->lastPt());
  } else {
    muf = dipole()->realEmissionME()->factorizationScale();
  }
  muf *= dipole()->realEmissionME()->factorizationScaleFactor();
  double pdfweight = 1.;
  if ( realCXComb()->mePartonData()[0]->coloured() &&
       dipole()->realEmissionME()->havePDFWeight1() )
    pdfweight *= dipole()->realEmissionME()->pdf1(muf,theExtrapolationX);
  if ( realCXComb()->mePartonData()[1]->coloured() &&
       dipole()->realEmissionME()->havePDFWeight2() )
    pdfweight *= dipole()->realEmissionME()->pdf2(muf,theExtrapolationX);
  return pdfweight;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ShowerApproximation::persistentOutput(PersistentOStream & os) const {
  os << theBornXComb << theRealXComb << theTildeXCombs << theDipole << theBelowCutoff
     << ounit(theFFPtCut,GeV) << ounit(theFIPtCut,GeV)
     << ounit(theIIPtCut,GeV) << theShowerScalesInSubtraction
     << theShowerScalesInSplitting
     << theRestrictPhasespace << theHardScaleFactor
     << theExtrapolationX;
}

void ShowerApproximation::persistentInput(PersistentIStream & is, int) {
  is >> theBornXComb >> theRealXComb >> theTildeXCombs >> theDipole >> theBelowCutoff
     >> iunit(theFFPtCut,GeV) >> iunit(theFIPtCut,GeV)
     >> iunit(theIIPtCut,GeV) >> theShowerScalesInSubtraction
     >> theShowerScalesInSplitting
     >> theRestrictPhasespace >> theHardScaleFactor
     >> theExtrapolationX;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<ShowerApproximation,HandlerBase>
  describeHerwigShowerApproximation("Herwig::ShowerApproximation", "HwMatchbox.so");

void ShowerApproximation::Init() {

  static ClassDocumentation<ShowerApproximation> documentation
    ("ShowerApproximation describes the shower emission to be used "
     "in NLO matching.");

  static Parameter<ShowerApproximation,Energy> interfaceFFPtCut
    ("FFPtCut",
     "Set the pt infrared cutoff",
     &ShowerApproximation::theFFPtCut, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximation,Energy> interfaceFIPtCut
    ("FIPtCut",
     "Set the pt infrared cutoff",
     &ShowerApproximation::theFIPtCut, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximation,Energy> interfaceIIPtCut
    ("IIPtCut",
     "Set the pt infrared cutoff",
     &ShowerApproximation::theIIPtCut, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Switch<ShowerApproximation,bool> interfaceShowerScalesInSubtraction
    ("ShowerScalesInSubtraction",
     "Switch on or off shower scales in the matching subtraction",
     &ShowerApproximation::theShowerScalesInSubtraction, true, false, false);
  static SwitchOption interfaceShowerScalesInSubtractionOn
    (interfaceShowerScalesInSubtraction,
     "On",
     "Use shower scales in the matching subtraction",
     true);
  static SwitchOption interfaceShowerScalesInSubtractionOff
    (interfaceShowerScalesInSubtraction,
     "Off",
     "Use hard process scales in the matching subtraction",
     false);

  static Switch<ShowerApproximation,bool> interfaceShowerScalesInSplitting
    ("ShowerScalesInSplitting",
     "Switch on or off shower scales in the splitting generation",
     &ShowerApproximation::theShowerScalesInSplitting, true, false, false);
  static SwitchOption interfaceShowerScalesInSplittingOn
    (interfaceShowerScalesInSplitting,
     "On",
     "Use shower scales in the matching splitting generation",
     true);
  static SwitchOption interfaceShowerScalesInSplittingOff
    (interfaceShowerScalesInSplitting,
     "Off",
     "Use hard process scales in the matching splitting generation",
     false);

  static Switch<ShowerApproximation,bool> interfaceRestrictPhasespace
    ("RestrictPhasespace",
     "Switch on or off phasespace restrictions",
     &ShowerApproximation::theRestrictPhasespace, true, false, false);
  static SwitchOption interfaceRestrictPhasespaceOn
    (interfaceRestrictPhasespace,
     "On",
     "Perform phasespace restrictions",
     true);
  static SwitchOption interfaceRestrictPhasespaceOff
    (interfaceRestrictPhasespace,
     "Off",
     "Do not perform phasespace restrictions",
     false);

  static Parameter<ShowerApproximation,double> interfaceHardScaleFactor
    ("HardScaleFactor",
     "The hard scale factor.",
     &ShowerApproximation::theHardScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximation,double> interfaceExtrapolationX
    ("ExtrapolationX",
     "The x from which on extrapolation should be performed.",
     &ShowerApproximation::theExtrapolationX, 0.65, 0.0, 1.0,
     false, false, Interface::limited);

}

