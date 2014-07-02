// -*- C++ -*-
//
// DipoleSplittingKernel.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleSplittingKernel class.
//

#include "DipoleSplittingKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DipoleSplittingKernel::DipoleSplittingKernel() 
  : HandlerBase(), theScreeningScale(0.0*GeV), 
    thePresamplingPoints(50000), theMaxtry(100000),
    theFreezeGrid(500000),
    theStrictLargeN(false), 
    theFactorizationScaleFactor(1.0),
    theRenormalizationScaleFactor(1.0),
    theRenormalizationScaleFreeze(1.*GeV), 
    theFactorizationScaleFreeze(1.*GeV) {}

DipoleSplittingKernel::~DipoleSplittingKernel() {}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleSplittingKernel::persistentOutput(PersistentOStream & os) const {
  os << theAlphaS << ounit(theScreeningScale,GeV) << theSplittingKinematics << thePDFRatio
     << thePresamplingPoints << theMaxtry << theFreezeGrid
     << theFlavour << theMCCheck << theStrictLargeN
     << theFactorizationScaleFactor
     << theRenormalizationScaleFactor
     << ounit(theRenormalizationScaleFreeze,GeV)
     << ounit(theFactorizationScaleFreeze,GeV);
}

void DipoleSplittingKernel::persistentInput(PersistentIStream & is, int) {
  is >> theAlphaS >> iunit(theScreeningScale,GeV) >> theSplittingKinematics >> thePDFRatio
     >> thePresamplingPoints >> theMaxtry >> theFreezeGrid
     >> theFlavour >> theMCCheck >> theStrictLargeN
     >> theFactorizationScaleFactor
     >> theRenormalizationScaleFactor
     >> iunit(theRenormalizationScaleFreeze,GeV)
     >> iunit(theFactorizationScaleFreeze,GeV);
}

double DipoleSplittingKernel::alphaPDF(const DipoleSplittingInfo& split) const {

  Energy pt = split.lastPt();

  Energy2 scale = sqr(pt) + sqr(theScreeningScale);

  Energy2 rScale = sqr(theRenormalizationScaleFactor)*scale;
  rScale = rScale > sqr(renormalizationScaleFreeze()) ? rScale : sqr(renormalizationScaleFreeze());

  Energy2 fScale = sqr(theFactorizationScaleFactor)*scale;
  fScale = fScale > sqr(factorizationScaleFreeze()) ? fScale : sqr(factorizationScaleFreeze());

  double ret = alphaS()->value(rScale) / (2.*Constants::pi);

  if ( split.index().initialStateEmitter() ) {
    assert(pdfRatio());
    ret *= 
      split.lastEmitterZ() * 
      (*pdfRatio())(split.index().emitterPDF(), fScale,
		    split.index().emitterData(),split.emitterData(),
		    split.emitterX(),split.lastEmitterZ());
  }

  if ( split.index().initialStateSpectator() ) {
    assert(pdfRatio());
    ret *= 
      split.lastSpectatorZ() * 
      (*pdfRatio())(split.index().spectatorPDF(), fScale,
		    split.index().spectatorData(),split.spectatorData(),
		    split.spectatorX(),split.lastSpectatorZ());
  }


  if ( ret < 0. )
    ret = 0.;

  return ret;

}

AbstractClassDescription<DipoleSplittingKernel> DipoleSplittingKernel::initDipoleSplittingKernel;
// Definition of the static class description member.

void DipoleSplittingKernel::Init() {

  static ClassDocumentation<DipoleSplittingKernel> documentation
    ("DipoleSplittingKernel is the base class for all kernels "
     "used within the dipole shower.");

  static Reference<DipoleSplittingKernel,AlphaSBase> interfaceAlphaS
    ("AlphaS",
     "The strong coupling to be used by this splitting kernel.",
     &DipoleSplittingKernel::theAlphaS, false, false, true, true, false);


  static Parameter<DipoleSplittingKernel,Energy> interfaceScreeningScale
    ("ScreeningScale",
     "A colour screening scale",
     &DipoleSplittingKernel::theScreeningScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


  static Reference<DipoleSplittingKernel,DipoleSplittingKinematics> interfaceSplittingKinematics
    ("SplittingKinematics",
     "The splitting kinematics to be used by this splitting kernel.",
     &DipoleSplittingKernel::theSplittingKinematics, false, false, true, false, false);


  static Reference<DipoleSplittingKernel,PDFRatio> interfacePDFRatio
    ("PDFRatio",
     "Set the optional PDF ratio object to evaluate this kernel",
     &DipoleSplittingKernel::thePDFRatio, false, false, true, true, false);

  static Parameter<DipoleSplittingKernel,unsigned long> interfacePresamplingPoints
    ("PresamplingPoints",
     "The number of points used to presample this kernel.",
     &DipoleSplittingKernel::thePresamplingPoints, 50000, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<DipoleSplittingKernel,unsigned long> interfaceMaxtry
    ("Maxtry",
     "The maximum number of attempts to generate a splitting.",
     &DipoleSplittingKernel::theMaxtry, 10000, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<DipoleSplittingKernel,unsigned long> interfaceFreezeGrid
    ("FreezeGrid",
     "",
     &DipoleSplittingKernel::theFreezeGrid, 500000, 1, 0,
     false, false, Interface::lowerlim);

  static Reference<DipoleSplittingKernel,ParticleData> interfaceFlavour
    ("Flavour",
     "Set the flavour to be produced if ambiguous.",
     &DipoleSplittingKernel::theFlavour, false, false, true, true, false);

  static Reference<DipoleSplittingKernel,DipoleMCCheck> interfaceMCCheck
    ("MCCheck",
     "[debug option] MCCheck",
     &DipoleSplittingKernel::theMCCheck, false, false, true, true, false);

  interfaceMCCheck.rank(-1);

  static Switch<DipoleSplittingKernel,bool> interfaceStrictLargeN
    ("StrictLargeN",
     "Work in a strict large-N limit.",
     &DipoleSplittingKernel::theStrictLargeN, false, false, false);
  static SwitchOption interfaceStrictLargeNOn
    (interfaceStrictLargeN,
     "On",
     "Replace C_F -> C_A/2 where present",
     true);
  static SwitchOption interfaceStrictLargeNOff
    (interfaceStrictLargeN,
     "Off",
     "Keep C_F=4/3",
     false);

  interfaceStrictLargeN.rank(-2);

  static Parameter<DipoleSplittingKernel,double> interfaceFactorizationScaleFactor
    ("FactorizationScaleFactor",
     "The factorization scale factor.",
     &DipoleSplittingKernel::theFactorizationScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  interfaceFactorizationScaleFactor.rank(-2);

  static Parameter<DipoleSplittingKernel,double> interfaceRenormalizationScaleFactor
    ("RenormalizationScaleFactor",
     "The renormalization scale factor.",
     &DipoleSplittingKernel::theRenormalizationScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  interfaceRenormalizationScaleFactor.rank(-2);

  static Parameter<DipoleSplittingKernel,Energy> interfaceRenormalizationScaleFreeze
    ("RenormalizationScaleFreeze",
     "The freezing scale for the renormalization scale.",
     &DipoleSplittingKernel::theRenormalizationScaleFreeze, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<DipoleSplittingKernel,Energy> interfaceFactorizationScaleFreeze
    ("FactorizationScaleFreeze",
     "The freezing scale for the factorization scale.",
     &DipoleSplittingKernel::theFactorizationScaleFreeze, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


}

