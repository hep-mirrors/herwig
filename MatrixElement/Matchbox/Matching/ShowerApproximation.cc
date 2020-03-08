// -*- C++ -*-
//
// ShowerApproximation.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerApproximation class.
//

#include "ShowerApproximation.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/TildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/InvertedTildeKinematics.h"

using namespace Herwig;

ShowerApproximation::ShowerApproximation() 
  : HandlerBase(),
    theExtrapolationX(1.0), theBelowCutoff(false),
    theFFPtCut(1.0*GeV), theFFScreeningScale(ZERO),
    theFIPtCut(1.0*GeV), theFIScreeningScale(ZERO),
    theIIPtCut(1.0*GeV), theIIScreeningScale(ZERO),
    theSafeCut(0.0*GeV),
    theRestrictPhasespace(true), theHardScaleFactor(1.0),
    theRenormalizationScaleFactor(1.0), theFactorizationScaleFactor(1.0),
    theRealEmissionScaleInSubtraction(showerScale), 
    theBornScaleInSubtraction(showerScale), 
    theEmissionScaleInSubtraction(showerScale), 
    theRealEmissionScaleInSplitting(showerScale), 
    theBornScaleInSplitting(showerScale), 
    theEmissionScaleInSplitting(showerScale),
    theRenormalizationScaleFreeze(1.*GeV),
    theFactorizationScaleFreeze(1.*GeV),
  maxPtIsMuF(false), theOpenZ(true) {}

ShowerApproximation::~ShowerApproximation() {}

void ShowerApproximation::setLargeNBasis() {
  assert(dipole()->realEmissionME()->matchboxAmplitude());
  if ( !dipole()->realEmissionME()->matchboxAmplitude()->treeAmplitudes() )
    return;
  if ( !theLargeNBasis ) {
    if ( !dipole()->realEmissionME()->matchboxAmplitude()->colourBasis() )
      throw Exception() << "ShowerApproximation::setLargeNBasis(): Expecting a colour basis object."
			<< Exception::runerror;
    theLargeNBasis = 
      dipole()->realEmissionME()->matchboxAmplitude()->colourBasis()->cloneMe();
    theLargeNBasis->clear();
    theLargeNBasis->doLargeN();
  }
}

void ShowerApproximation::setDipole(Ptr<SubtractionDipole>::tptr dip) { 
  theDipole = dip;
  setLargeNBasis();
}

Ptr<SubtractionDipole>::tptr ShowerApproximation::dipole() const { return theDipole; }

Ptr<TildeKinematics>::tptr
ShowerApproximation::showerTildeKinematics() const {
  return Ptr<TildeKinematics>::tptr();
}

Ptr<InvertedTildeKinematics>::tptr 
ShowerApproximation::showerInvertedTildeKinematics() const {
  return Ptr<InvertedTildeKinematics>::tptr();
}

void ShowerApproximation::checkCutoff() {
  assert(!showerTildeKinematics());
}

void ShowerApproximation::getShowerVariables() {

  // check for the cutoff
  dipole()->isAboveCutoff(isAboveCutoff());

  // get the hard scale
  dipole()->showerHardScale(hardScale());

  // set the shower scale and variables for completeness
  dipole()->showerScale(dipole()->lastPt());
  dipole()->showerParameters().resize(1);
  dipole()->showerParameters()[0] = dipole()->lastZ();

  // check for phase space
  dipole()->isInShowerPhasespace(isInShowerPhasespace());

}

bool ShowerApproximation::isAboveCutoff() const {

  if ( dipole()->bornEmitter() > 1 &&
       dipole()->bornSpectator() > 1 ) {
    return dipole()->lastPt() >= max(ffPtCut(),safeCut());
  } else if ( ( dipole()->bornEmitter() > 1 &&
		dipole()->bornSpectator() < 2 ) ||
	      ( dipole()->bornEmitter() < 2 &&
		dipole()->bornSpectator() > 1 ) ) {
    return dipole()->lastPt() >= max(fiPtCut(),safeCut());
  } else {
    assert(dipole()->bornEmitter() < 2 &&
	   dipole()->bornSpectator() < 2);
    return dipole()->lastPt() >= max(iiPtCut(),safeCut());
  }

  return true;

}

Energy ShowerApproximation::hardScale() const {
  if ( !maxPtIsMuF ) {
    if ( !bornCXComb()->mePartonData()[0]->coloured() &&
	 !bornCXComb()->mePartonData()[1]->coloured() ) {
      Energy maxPt = (bornCXComb()->meMomenta()[0] + bornCXComb()->meMomenta()[1]).m();
      maxPt *= hardScaleFactor();
      return maxPt;
    }
    Energy maxPt = generator()->maximumCMEnergy();
    vector<Lorentz5Momentum>::const_iterator p = 
      bornCXComb()->meMomenta().begin() + 2;
    cPDVector::const_iterator pp = 
      bornCXComb()->mePartonData().begin() + 2;
    for ( ; p != bornCXComb()->meMomenta().end(); ++p, ++pp )
      if ( (**pp).coloured() )
	maxPt = min(maxPt,p->mt());
    if ( maxPt == generator()->maximumCMEnergy() )
      maxPt = (bornCXComb()->meMomenta()[0] + bornCXComb()->meMomenta()[1]).m();
    maxPt *= hardScaleFactor();
    return maxPt;
  } else {
    return hardScaleFactor()*sqrt(bornCXComb()->lastShowerScale());
  }
}

bool ShowerApproximation::isInShowerPhasespace() const {

  if ( !dipole()->isAboveCutoff() )
    return false;
  if ( !restrictPhasespace() )
    return true;

  InvertedTildeKinematics& kinematics =
    const_cast<InvertedTildeKinematics&>(*dipole()->invertedTildeKinematics());
  tcStdXCombPtr tmpreal = kinematics.realXComb();
  tcStdXCombPtr tmpborn = kinematics.bornXComb();
  Ptr<SubtractionDipole>::tptr tmpdip = kinematics.dipole();

  Energy hard = dipole()->showerHardScale();
  Energy pt = dipole()->lastPt();
  double z = dipole()->lastZ();

  pair<double,double> zbounds(0.,1.);

  kinematics.dipole(const_ptr_cast<Ptr<SubtractionDipole>::tptr>(theDipole));
  kinematics.prepare(realCXComb(),bornCXComb());

  if ( pt > hard ) {
    kinematics.dipole(tmpdip);
    kinematics.prepare(tmpreal,tmpborn);
    return false;
  }

  try {
    zbounds = kinematics.zBounds(pt,openZ() ? kinematics.ptMax() : hard);
  } catch(...) {
    kinematics.dipole(tmpdip);
    kinematics.prepare(tmpreal,tmpborn);
    throw;
  }
  kinematics.dipole(tmpdip);
  kinematics.prepare(tmpreal,tmpborn);

  return z > zbounds.first && z < zbounds.second;

}

Energy2 ShowerApproximation::showerEmissionScale() const {

  Energy2 mur = sqr(dipole()->lastPt());

  if ( dipole()->bornEmitter() > 1 &&
       dipole()->bornSpectator() > 1 ) {
    return mur + sqr(ffScreeningScale());
  } else if ( ( dipole()->bornEmitter() > 1 &&
		dipole()->bornSpectator() < 2 ) ||
	      ( dipole()->bornEmitter() < 2 &&
		dipole()->bornSpectator() > 1 ) ) {
    return mur + sqr(fiScreeningScale());
  } else {
    assert(dipole()->bornEmitter() < 2 &&
	   dipole()->bornSpectator() < 2);
    return mur + sqr(iiScreeningScale());
  }

  return mur;

}

Energy2 ShowerApproximation::bornRenormalizationScale() const {
  return 
    sqr(dipole()->underlyingBornME()->renormalizationScaleFactor()) *
    dipole()->underlyingBornME()->renormalizationScale();
}

Energy2 ShowerApproximation::bornFactorizationScale() const {
  return 
    sqr(dipole()->underlyingBornME()->factorizationScaleFactor()) *
    dipole()->underlyingBornME()->factorizationScale();
}

Energy2 ShowerApproximation::realRenormalizationScale() const {
  return 
    sqr(dipole()->realEmissionME()->renormalizationScaleFactor()) *
    dipole()->realEmissionME()->renormalizationScale();
}

Energy2 ShowerApproximation::realFactorizationScale() const {
  return 
    sqr(dipole()->realEmissionME()->factorizationScaleFactor()) *
    dipole()->realEmissionME()->factorizationScale();
}

double ShowerApproximation::bornPDFWeight(Energy2 muf) const {
  if ( !bornCXComb()->mePartonData()[0]->coloured() &&
       !bornCXComb()->mePartonData()[1]->coloured() )
    return 1.;
  if ( muf < sqr(theFactorizationScaleFreeze) )
    muf = sqr(theFactorizationScaleFreeze);
  double pdfweight = 1.;
  if ( bornCXComb()->mePartonData()[0]->coloured() &&
       dipole()->underlyingBornME()->havePDFWeight1() )
    pdfweight *= dipole()->underlyingBornME()->pdf1(muf,theExtrapolationX);
  if ( bornCXComb()->mePartonData()[1]->coloured() &&
       dipole()->underlyingBornME()->havePDFWeight2() )
    pdfweight *= dipole()->underlyingBornME()->pdf2(muf,theExtrapolationX);
  return pdfweight;
}

double ShowerApproximation::realPDFWeight(Energy2 muf) const {
  if ( !realCXComb()->mePartonData()[0]->coloured() &&
       !realCXComb()->mePartonData()[1]->coloured() )
    return 1.;
  if ( muf < sqr(theFactorizationScaleFreeze) )
    muf = sqr(theFactorizationScaleFreeze);
  double pdfweight = 1.;
  if ( realCXComb()->mePartonData()[0]->coloured() &&
       dipole()->realEmissionME()->havePDFWeight1() )
    pdfweight *= dipole()->realEmissionME()->pdf1(muf,theExtrapolationX);
  if ( realCXComb()->mePartonData()[1]->coloured() &&
       dipole()->realEmissionME()->havePDFWeight2() )
    pdfweight *= dipole()->realEmissionME()->pdf2(muf,theExtrapolationX);
  return pdfweight;
}

double ShowerApproximation::scaleWeight(int rScale, int bScale, int eScale) const {

  double emissionAlpha = 1.;
  Energy2 emissionScale = ZERO;
  Energy2 showerscale = ZERO;
  if ( eScale == showerScale || bScale == showerScale || eScale == showerScale ) {
    showerscale = showerRenormalizationScale();
    if ( showerscale < sqr(theRenormalizationScaleFreeze) )
      showerscale = sqr(theFactorizationScaleFreeze);
  }
  if ( eScale == showerScale ) {
    emissionAlpha = SM().alphaS(showerscale);
    emissionScale = showerFactorizationScale();
  } else if ( eScale == realScale ) {
    emissionAlpha = dipole()->realEmissionME()->lastXComb().lastAlphaS();
    emissionScale = dipole()->realEmissionME()->lastScale();
  } else if ( eScale == bornScale ) {
    emissionAlpha = dipole()->underlyingBornME()->lastXComb().lastAlphaS();
    emissionScale = dipole()->underlyingBornME()->lastScale();
  }
  double emissionPDF = realPDFWeight(emissionScale);

  double couplingFactor = 1.;
  if ( bScale != rScale ) {
    double bornAlpha = 1.;
    if ( bScale == showerScale ) {
      bornAlpha = SM().alphaS(showerscale);
    } else if ( bScale == realScale ) {
      bornAlpha = dipole()->realEmissionME()->lastXComb().lastAlphaS();
    } else if ( bScale == bornScale ) {
      bornAlpha = dipole()->underlyingBornME()->lastXComb().lastAlphaS();
    }
    double realAlpha = 1.;
    if ( rScale == showerScale ) {
      realAlpha = SM().alphaS(showerscale);
    } else if ( rScale == realScale ) {
      realAlpha = dipole()->realEmissionME()->lastXComb().lastAlphaS();
    } else if ( rScale == bornScale ) {
      realAlpha = dipole()->underlyingBornME()->lastXComb().lastAlphaS();
    }
    couplingFactor *=
      pow(realAlpha/bornAlpha,(double)(dipole()->underlyingBornME()->orderInAlphaS()));
  }

  Energy2 hardScale = ZERO;
  if ( bScale == showerScale ) {
    hardScale = showerFactorizationScale();
  } else if ( bScale == realScale ) {
    hardScale = dipole()->realEmissionME()->lastScale();
  } else if ( bScale == bornScale ) {
    hardScale = dipole()->underlyingBornME()->lastScale();
  }
  double bornPDF = bornPDFWeight(hardScale);

  if ( abs(bornPDF) < 1e-8 )
    bornPDF = 0.0;

  if ( abs(emissionPDF) < 1e-8 )
    emissionPDF = 0.0;

  if ( emissionPDF == 0.0 || bornPDF == 0.0 )
    return 0.0;

  double pdfRatio = emissionPDF/bornPDF;
  pdfRatio = min(abs(pdfRatio),100000.);

  return
    emissionAlpha * pdfRatio * couplingFactor;
    
}

double ShowerApproximation::channelWeight(int emitter, int emission, 
					  int spectator, int) const {
  double cfac = 1.;
  double Nc = generator()->standardModel()->Nc();
  if (realCXComb()->mePartonData()[emitter]->iColour() == PDT::Colour8){
    if (realCXComb()->mePartonData()[emission]->iColour() == PDT::Colour8)
      cfac = Nc;
    else if ( realCXComb()->mePartonData()[emission]->iColour() == PDT::Colour3 ||
              realCXComb()->mePartonData()[emission]->iColour() == PDT::Colour3bar)
      cfac = 0.5;
    else assert(false);
  }
  else if ((realCXComb()->mePartonData()[emitter] ->iColour() == PDT::Colour3 ||
            realCXComb()->mePartonData()[emitter] ->iColour() == PDT::Colour3bar))
    cfac = (sqr(Nc)-1.)/(2.*Nc);
  else assert(false);
  // do the most simple thing for the time being; needs fixing later
  if ( realCXComb()->mePartonData()[emission]->id() == ParticleID::g ) {
    Energy2 pipk = 
      realCXComb()->meMomenta()[emitter] * realCXComb()->meMomenta()[spectator];
    Energy2 pipj = 
      realCXComb()->meMomenta()[emitter] * realCXComb()->meMomenta()[emission];
    Energy2 pjpk = 
      realCXComb()->meMomenta()[emission] * realCXComb()->meMomenta()[spectator];
    return cfac *GeV2 * pipk / ( pipj * ( pipj + pjpk ) );
  }
  return
    cfac * GeV2 / (realCXComb()->meMomenta()[emitter] * realCXComb()->meMomenta()[emission]);
}

double ShowerApproximation::channelWeight() const {
  double currentChannel = channelWeight(dipole()->realEmitter(),
					dipole()->realEmission(),
					dipole()->realSpectator(),
					dipole()->bornEmitter());
  if ( currentChannel == 0. )
    return 0.;
  double sum = 0.;
  for ( vector<Ptr<SubtractionDipole>::tptr>::const_iterator dip =
	  dipole()->partnerDipoles().begin();
	dip != dipole()->partnerDipoles().end(); ++dip )
    sum += channelWeight((**dip).realEmitter(),
			 (**dip).realEmission(),
			 (**dip).realSpectator(),
			 (**dip).bornEmitter());
  assert(sum > 0.0);
  return currentChannel / sum;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void ShowerApproximation::doinit() {
  if ( profileScales() ) {
    if ( profileScales()->unrestrictedPhasespace() &&
	 restrictPhasespace() ) {
      generator()->log()
	<< "ShowerApproximation warning: The scale profile chosen requires an unrestricted phase space,\n"
	<< "however, the phase space was set to be restricted. Will switch to unrestricted phase space.\n"
	<< flush;
      restrictPhasespace(false);
    }
  }
  HandlerBase::doinit();
}

void ShowerApproximation::persistentOutput(PersistentOStream & os) const {
  os << theLargeNBasis
     << theBornXComb << theRealXComb << theTildeXCombs << theDipole << theBelowCutoff
     << ounit(theFFPtCut,GeV) << ounit(theFFScreeningScale,GeV) 
     << ounit(theFIPtCut,GeV) << ounit(theFIScreeningScale,GeV) 
     << ounit(theIIPtCut,GeV) << ounit(theIIScreeningScale,GeV)
     << ounit(theSafeCut,GeV) 
     << theRestrictPhasespace << theHardScaleFactor
     << theRenormalizationScaleFactor << theFactorizationScaleFactor
     << theExtrapolationX
     << theRealEmissionScaleInSubtraction << theBornScaleInSubtraction
     << theEmissionScaleInSubtraction << theRealEmissionScaleInSplitting
     << theBornScaleInSplitting << theEmissionScaleInSplitting
     << ounit(theRenormalizationScaleFreeze,GeV)
     << ounit(theFactorizationScaleFreeze,GeV) << maxPtIsMuF
     << theHardScaleProfile << theOpenZ;
}

void ShowerApproximation::persistentInput(PersistentIStream & is, int) {
  is >> theLargeNBasis
     >> theBornXComb >> theRealXComb >> theTildeXCombs >> theDipole >> theBelowCutoff
     >> iunit(theFFPtCut,GeV) >> iunit(theFFScreeningScale,GeV) 
     >> iunit(theFIPtCut,GeV) >> iunit(theFIScreeningScale,GeV) 
     >> iunit(theIIPtCut,GeV) >> iunit(theIIScreeningScale,GeV) 
     >> iunit(theSafeCut,GeV)
     >> theRestrictPhasespace >> theHardScaleFactor
     >> theRenormalizationScaleFactor >> theFactorizationScaleFactor
     >> theExtrapolationX
     >> theRealEmissionScaleInSubtraction >> theBornScaleInSubtraction
     >> theEmissionScaleInSubtraction >> theRealEmissionScaleInSplitting
     >> theBornScaleInSplitting >> theEmissionScaleInSplitting
     >> iunit(theRenormalizationScaleFreeze,GeV)
     >> iunit(theFactorizationScaleFreeze,GeV) >> maxPtIsMuF
     >> theHardScaleProfile >> theOpenZ;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<ShowerApproximation,HandlerBase>
  describeHerwigShowerApproximation("Herwig::ShowerApproximation", "Herwig.so");

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
    
  static Parameter<ShowerApproximation,Energy> interfaceSafeCut
    ("SafeCut",
     "Set the enhanced infrared cutoff for the Matching.",
     &ShowerApproximation::theSafeCut, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximation,Energy> interfaceFFScreeningScale
    ("FFScreeningScale",
     "Set the screening scale",
     &ShowerApproximation::theFFScreeningScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximation,Energy> interfaceFIScreeningScale
    ("FIScreeningScale",
     "Set the screening scale",
     &ShowerApproximation::theFIScreeningScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximation,Energy> interfaceIIScreeningScale
    ("IIScreeningScale",
     "Set the screening scale",
     &ShowerApproximation::theIIScreeningScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Switch<ShowerApproximation,bool> interfaceRestrictPhasespace
    ("RestrictPhasespace",
     "Switch on or off phasespace restrictions",
     &ShowerApproximation::theRestrictPhasespace, true, false, false);
  static SwitchOption interfaceRestrictPhasespaceYes
    (interfaceRestrictPhasespace,
     "Yes",
     "Perform phasespace restrictions",
     true);
  static SwitchOption interfaceRestrictPhasespaceNo
    (interfaceRestrictPhasespace,
     "No",
     "Do not perform phasespace restrictions",
     false);

  static Parameter<ShowerApproximation,double> interfaceHardScaleFactor
    ("HardScaleFactor",
     "The hard scale factor.",
     &ShowerApproximation::theHardScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximation,double> interfaceRenormalizationScaleFactor
    ("RenormalizationScaleFactor",
     "The hard scale factor.",
     &ShowerApproximation::theRenormalizationScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximation,double> interfaceFactorizationScaleFactor
    ("FactorizationScaleFactor",
     "The hard scale factor.",
     &ShowerApproximation::theFactorizationScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximation,double> interfaceExtrapolationX
    ("ExtrapolationX",
     "The x from which on extrapolation should be performed.",
     &ShowerApproximation::theExtrapolationX, 1.0, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<ShowerApproximation,int> interfaceRealEmissionScaleInSubtraction
    ("RealEmissionScaleInSubtraction",
     "Set the scale choice for the real emission cross section in the matching subtraction.",
     &ShowerApproximation::theRealEmissionScaleInSubtraction, showerScale, false, false);
  static SwitchOption interfaceRealEmissionScaleInSubtractionRealScale
    (interfaceRealEmissionScaleInSubtraction,
     "RealScale",
     "Use the real emission scale.",
     realScale);
  static SwitchOption interfaceRealEmissionScaleInSubtractionBornScale
    (interfaceRealEmissionScaleInSubtraction,
     "BornScale",
     "Use the Born scale.",
     bornScale);
  static SwitchOption interfaceRealEmissionScaleInSubtractionShowerScale
    (interfaceRealEmissionScaleInSubtraction,
     "ShowerScale",
     "Use the shower scale",
     showerScale);

  interfaceRealEmissionScaleInSubtraction.rank(-1);

  static Switch<ShowerApproximation,int> interfaceBornScaleInSubtraction
    ("BornScaleInSubtraction",
     "Set the scale choice for the Born cross section in the matching subtraction.",
     &ShowerApproximation::theBornScaleInSubtraction, showerScale, false, false);
  static SwitchOption interfaceBornScaleInSubtractionRealScale
    (interfaceBornScaleInSubtraction,
     "RealScale",
     "Use the real emission scale.",
     realScale);
  static SwitchOption interfaceBornScaleInSubtractionBornScale
    (interfaceBornScaleInSubtraction,
     "BornScale",
     "Use the Born scale.",
     bornScale);
  static SwitchOption interfaceBornScaleInSubtractionShowerScale
    (interfaceBornScaleInSubtraction,
     "ShowerScale",
     "Use the shower scale",
     showerScale);

  interfaceBornScaleInSubtraction.rank(-1);

  static Switch<ShowerApproximation,int> interfaceEmissionScaleInSubtraction
    ("EmissionScaleInSubtraction",
     "Set the scale choice for the emission in the matching subtraction.",
     &ShowerApproximation::theEmissionScaleInSubtraction, showerScale, false, false);
  static SwitchOption interfaceEmissionScaleInSubtractionRealScale
    (interfaceEmissionScaleInSubtraction,
     "RealScale",
     "Use the real emission scale.",
     realScale);
  static SwitchOption interfaceEmissionScaleInSubtractionEmissionScale
    (interfaceEmissionScaleInSubtraction,
     "BornScale",
     "Use the Born scale.",
     bornScale);
  static SwitchOption interfaceEmissionScaleInSubtractionShowerScale
    (interfaceEmissionScaleInSubtraction,
     "ShowerScale",
     "Use the shower scale",
     showerScale);

  interfaceEmissionScaleInSubtraction.rank(-1);

  static Switch<ShowerApproximation,int> interfaceRealEmissionScaleInSplitting
    ("RealEmissionScaleInSplitting",
     "Set the scale choice for the real emission cross section in the splitting.",
     &ShowerApproximation::theRealEmissionScaleInSplitting, showerScale, false, false);
  static SwitchOption interfaceRealEmissionScaleInSplittingRealScale
    (interfaceRealEmissionScaleInSplitting,
     "RealScale",
     "Use the real emission scale.",
     realScale);
  static SwitchOption interfaceRealEmissionScaleInSplittingBornScale
    (interfaceRealEmissionScaleInSplitting,
     "BornScale",
     "Use the Born scale.",
     bornScale);
  static SwitchOption interfaceRealEmissionScaleInSplittingShowerScale
    (interfaceRealEmissionScaleInSplitting,
     "ShowerScale",
     "Use the shower scale",
     showerScale);

  interfaceRealEmissionScaleInSplitting.rank(-1);

  static Switch<ShowerApproximation,int> interfaceBornScaleInSplitting
    ("BornScaleInSplitting",
     "Set the scale choice for the Born cross section in the splitting.",
     &ShowerApproximation::theBornScaleInSplitting, showerScale, false, false);
  static SwitchOption interfaceBornScaleInSplittingRealScale
    (interfaceBornScaleInSplitting,
     "RealScale",
     "Use the real emission scale.",
     realScale);
  static SwitchOption interfaceBornScaleInSplittingBornScale
    (interfaceBornScaleInSplitting,
     "BornScale",
     "Use the Born scale.",
     bornScale);
  static SwitchOption interfaceBornScaleInSplittingShowerScale
    (interfaceBornScaleInSplitting,
     "ShowerScale",
     "Use the shower scale",
     showerScale);

  interfaceBornScaleInSplitting.rank(-1);

  static Switch<ShowerApproximation,int> interfaceEmissionScaleInSplitting
    ("EmissionScaleInSplitting",
     "Set the scale choice for the emission in the splitting.",
     &ShowerApproximation::theEmissionScaleInSplitting, showerScale, false, false);
  static SwitchOption interfaceEmissionScaleInSplittingRealScale
    (interfaceEmissionScaleInSplitting,
     "RealScale",
     "Use the real emission scale.",
     realScale);
  static SwitchOption interfaceEmissionScaleInSplittingEmissionScale
    (interfaceEmissionScaleInSplitting,
     "BornScale",
     "Use the Born scale.",
     bornScale);
  static SwitchOption interfaceEmissionScaleInSplittingShowerScale
    (interfaceEmissionScaleInSplitting,
     "ShowerScale",
     "Use the shower scale",
     showerScale);

  interfaceEmissionScaleInSplitting.rank(-1);

  static Parameter<ShowerApproximation,Energy> interfaceRenormalizationScaleFreeze
    ("RenormalizationScaleFreeze",
     "The freezing scale for the renormalization scale.",
     &ShowerApproximation::theRenormalizationScaleFreeze, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  interfaceRenormalizationScaleFreeze.rank(-1);

  static Parameter<ShowerApproximation,Energy> interfaceFactorizationScaleFreeze
    ("FactorizationScaleFreeze",
     "The freezing scale for the factorization scale.",
     &ShowerApproximation::theFactorizationScaleFreeze, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  interfaceFactorizationScaleFreeze.rank(-1);

  static Reference<ShowerApproximation,HardScaleProfile> interfaceHardScaleProfile
    ("HardScaleProfile",
     "The hard scale profile to use.",
     &ShowerApproximation::theHardScaleProfile, false, false, true, true, false);

  static Reference<ShowerApproximation,ColourBasis> interfaceLargeNBasis
    ("LargeNBasis",
     "Set the large-N colour basis implementation.",
     &ShowerApproximation::theLargeNBasis, false, false, true, true, false);

  interfaceLargeNBasis.rank(-1);

  static Switch<ShowerApproximation,bool> interfaceMaxPtIsMuF
    ("MaxPtIsMuF",
     "",
     &ShowerApproximation::maxPtIsMuF, false, false, false);
  static SwitchOption interfaceMaxPtIsMuFYes
    (interfaceMaxPtIsMuF,
     "Yes",
     "",
     true);
  static SwitchOption interfaceMaxPtIsMuFNo
    (interfaceMaxPtIsMuF,
     "No",
     "",
     false);

}

