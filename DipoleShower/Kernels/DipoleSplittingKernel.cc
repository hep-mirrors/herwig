// -*- C++ -*-
//
// DipoleSplittingKernel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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

#include "Herwig/Shower/ShowerHandler.h"

using namespace Herwig;

DipoleSplittingKernel::DipoleSplittingKernel() 
  : HandlerBase(), theScreeningScale(0.0*GeV), 
    thePresamplingPoints(2000), theMaxtry(100000),
    theFreezeGrid(500000),
    theDetuning(1.0),
    theStrictLargeN(false), 
    theFactorizationScaleFactor(1.0),
    theRenormalizationScaleFactor(1.0),
    theRenormalizationScaleFreeze(1.*GeV), 
    theFactorizationScaleFreeze(1.*GeV),
    theVirtualitySplittingScale(false),
    presampling(false) {}

DipoleSplittingKernel::~DipoleSplittingKernel() {}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleSplittingKernel::persistentOutput(PersistentOStream & os) const {
  os << theAlphaS << ounit(theScreeningScale,GeV) << theSplittingKinematics << thePDFRatio
     << thePresamplingPoints << theMaxtry << theFreezeGrid << theDetuning
     << theFlavour << theMCCheck << theStrictLargeN
     << theFactorizationScaleFactor
     << theRenormalizationScaleFactor
     << ounit(theRenormalizationScaleFreeze,GeV)
     << ounit(theFactorizationScaleFreeze,GeV)
     << theVirtualitySplittingScale;
}

void DipoleSplittingKernel::persistentInput(PersistentIStream & is, int) {
  is >> theAlphaS >> iunit(theScreeningScale,GeV) >> theSplittingKinematics >> thePDFRatio
     >> thePresamplingPoints >> theMaxtry >> theFreezeGrid >> theDetuning
     >> theFlavour >> theMCCheck >> theStrictLargeN
     >> theFactorizationScaleFactor
     >> theRenormalizationScaleFactor
     >> iunit(theRenormalizationScaleFreeze,GeV)
     >> iunit(theFactorizationScaleFreeze,GeV)
     >> theVirtualitySplittingScale;
}

double DipoleSplittingKernel::alphaPDF(const DipoleSplittingInfo& split,
				       Energy optScale,
				       double rScaleFactor,
				       double fScaleFactor) const {

  Energy pt = optScale == ZERO ? split.lastPt() : optScale;

  Energy2 scale = ZERO;
  if ( !virtualitySplittingScale() ) {
    scale = sqr(pt) + sqr(theScreeningScale);
  } else {
    scale = sqr(splittingKinematics()->QFromPt(pt,split)) + sqr(theScreeningScale);
  }

  Energy2 rScale = sqr(theRenormalizationScaleFactor*rScaleFactor)*scale;
  rScale = rScale > sqr(renormalizationScaleFreeze()) ? rScale : sqr(renormalizationScaleFreeze());

  Energy2 fScale = sqr(theFactorizationScaleFactor*fScaleFactor)*scale;
  fScale = fScale > sqr(factorizationScaleFreeze()) ? fScale : sqr(factorizationScaleFreeze());

  double alphas = 1.0;
  double pdf = 1.0;

  // check if we are potentially reweighting and cache evaluations
  bool evaluatePDF = true;
  bool evaluateAlphaS = true;
  bool variations = 
    !ShowerHandler::currentHandler()->showerVariations().empty() &&
    !presampling;
  if ( variations ) {
    /*
    cerr << "--------------------------------------------------------------------------------\n"
	 << flush;
    cerr << "variations ... central scale/GeV = "
	 <<  sqrt(scale/GeV2) << " r = "
	 << rScaleFactor << " f = " << fScaleFactor << "\n"
	 << flush;
    */
    map<double,double>::const_iterator pit = thePDFCache.find(fScaleFactor);
    evaluatePDF = (pit == thePDFCache.end());
    if ( !evaluatePDF ) {
      //cerr << "PDF is cached: ";
      pdf = pit->second;
      //cerr << pdf << "\n" << flush;
    }
    map<double,double>::const_iterator ait = theAlphaSCache.find(rScaleFactor);
    evaluateAlphaS = (ait == theAlphaSCache.end());
    if ( !evaluateAlphaS ) {
      //cerr << "alphas is cached: ";
      alphas = ait->second;
      //cerr << alphas << "\n" << flush;
    }
  }

  if ( evaluateAlphaS )
    alphas = alphaS()->value(rScale);

  if ( evaluatePDF ) {
    if ( split.index().initialStateEmitter() ) {
      assert(pdfRatio());
      pdf *= 
	split.lastEmitterZ() * 
	(*pdfRatio())(split.index().emitterPDF(), fScale,
		      split.index().emitterData(),split.emitterData(),
		      split.emitterX(),split.lastEmitterZ());
    }

    if ( split.index().initialStateSpectator() ) {
      assert(pdfRatio());
      pdf *= 
	split.lastSpectatorZ() * 
	(*pdfRatio())(split.index().spectatorPDF(), fScale,
		      split.index().spectatorData(),split.spectatorData(),
		      split.spectatorX(),split.lastSpectatorZ());
    }
  }

  if ( evaluatePDF && variations ) {
    //cerr << "caching PDF = " << pdf << "\n" << flush;
    thePDFCache[fScaleFactor] = pdf;
  }

  if ( evaluateAlphaS && variations ) {
    //cerr << "caching alphas = " << alphas << "\n" << flush;
    theAlphaSCache[rScaleFactor] = alphas;
  }

  double ret = alphas*pdf / (2.*Constants::pi);

  if ( ret < 0. )
    ret = 0.;

  return ret;

}

void DipoleSplittingKernel::accept(const DipoleSplittingInfo& split,
				   double, double,
				   map<string,double>& weights) const {
  if ( ShowerHandler::currentHandler()->showerVariations().empty() )
    return;
  double reference = alphaPDF(split);
  assert(reference > 0.);
  for ( map<string,ShowerHandler::ShowerVariation>::const_iterator var =
	  ShowerHandler::currentHandler()->showerVariations().begin();
	var != ShowerHandler::currentHandler()->showerVariations().end(); ++var ) {
    if ( ( ShowerHandler::currentHandler()->firstInteraction() && var->second.firstInteraction ) ||
	 ( !ShowerHandler::currentHandler()->firstInteraction() && var->second.secondaryInteractions ) ) {
      /*
      cerr << "reweighting in accept for: "
	   << var->first
	   << " using "
	   << var->second.renormalizationScaleFactor << " "
	   << var->second.factorizationScaleFactor << " "
	   << var->second.firstInteraction << " "
	   << var->second.secondaryInteractions << "\n" << flush;
      */
      double varied = alphaPDF(split,ZERO,
			       var->second.renormalizationScaleFactor,
			       var->second.factorizationScaleFactor);
      if ( varied != reference ) {
	map<string,double>::iterator wi = weights.find(var->first);
	if ( wi != weights.end() )
	  wi->second *= varied/reference;
	else
	  weights[var->first] = varied/reference;
      }
    }
  }
}

void DipoleSplittingKernel::veto(const DipoleSplittingInfo& split,
				 double p, double r,
				 map<string,double>& weights) const {
  if ( ShowerHandler::currentHandler()->showerVariations().empty() )
    return;
  double reference = alphaPDF(split);
  // this is dangerous, but we have no other choice currently -- need to
  // carefully check for the effects; the assumption is that if the central
  // one ius zero, then so will be the variations.
  if ( reference == 0.0 )
    return;
  for ( map<string,ShowerHandler::ShowerVariation>::const_iterator var =
	  ShowerHandler::currentHandler()->showerVariations().begin();
	var != ShowerHandler::currentHandler()->showerVariations().end(); ++var ) {
    if ( ( ShowerHandler::currentHandler()->firstInteraction() && var->second.firstInteraction ) ||
	 ( !ShowerHandler::currentHandler()->firstInteraction() && var->second.secondaryInteractions ) ) {
      /*
      cerr << "reweighting in veto for: "
	   << var->first
	   << " using "
	   << var->second.renormalizationScaleFactor << " "
	   << var->second.factorizationScaleFactor << " "
	   << var->second.firstInteraction << " "
	   << var->second.secondaryInteractions << "\n" << flush;
      */
      double varied = alphaPDF(split,ZERO,
			       var->second.renormalizationScaleFactor,
			       var->second.factorizationScaleFactor);
      if ( varied != reference ) {
	map<string,double>::iterator wi = weights.find(var->first);
	if ( wi != weights.end() )
	  wi->second *= (r - varied*p/reference) / (r-p);
	else
	  weights[var->first] = (r - varied*p/reference) / (r-p);
      }
    }
  }
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
     &DipoleSplittingKernel::thePresamplingPoints, 2000, 1, 0,
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

  static Switch<DipoleSplittingKernel,bool> interfaceVirtualitySplittingScale
    ("VirtualitySplittingScale",
     "Use the virtuality as the splitting scale.",
     &DipoleSplittingKernel::theVirtualitySplittingScale, false, false, false);
  static SwitchOption interfaceVirtualitySplittingScaleYes
    (interfaceVirtualitySplittingScale,
     "Yes",
     "Use vrituality.",
     true);
  static SwitchOption interfaceVirtualitySplittingScaleNo
    (interfaceVirtualitySplittingScale,
     "No",
     "Use transverse momentum.",
     false);

  static Parameter<DipoleSplittingKernel,double> interfaceDetuning
    ("Detuning",
     "A value to detune the overestimate kernel.",
     &DipoleSplittingKernel::theDetuning, 1.0, 1.0, 0,
     false, false, Interface::lowerlim);

}

