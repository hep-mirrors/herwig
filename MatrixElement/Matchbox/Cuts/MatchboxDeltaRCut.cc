// -*- C++ -*-
//
// MatchboxDeltaRCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxDeltaRCut class.
//

#include "MatchboxDeltaRCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxDeltaRCut::MatchboxDeltaRCut() 
  : theDeltaRMin(0.0), theDeltaRMax(Constants::MaxRapidity), 
    theDeltaYMin(0.0), theDeltaYMax(Constants::MaxRapidity),
    theDeltaPhiMin(0.0), theDeltaPhiMax(2.0*Constants::pi) {}

MatchboxDeltaRCut::~MatchboxDeltaRCut() {}

IBPtr MatchboxDeltaRCut::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxDeltaRCut::fullclone() const {
  return new_ptr(*this);
}

bool MatchboxDeltaRCut::passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
		                 LorentzMomentum pi, LorentzMomentum pj,
		                 bool inci, bool incj) const {

  bool match = false;
  if ( theFirstMatcher->check(*pitype) && theSecondMatcher->check(*pjtype) ) match = true;
  if ( theFirstMatcher->check(*pjtype) && theSecondMatcher->check(*pitype) ) match = true;
  if ( !match ||
       (theDeltaRMin == 0.0 && theDeltaRMax == Constants::MaxRapidity && 
	theDeltaYMin == 0.0 && theDeltaYMax == Constants::MaxRapidity &&
	theDeltaPhiMin == 0.0 && theDeltaPhiMax == 2.0*Constants::pi) ) return true;
  if ( inci || incj ) return true;

  double weight = 1.0;

  double dY = abs(pi.rapidity() - pj.rapidity());
  double dPhi = abs(pi.phi() - pj.phi());
  if ( dPhi > Constants::pi ) dPhi = 2.0*Constants::pi - dPhi;
  double dR = sqrt(sqr(dY) + sqr(dPhi));
  if ( !parent->isInside<CutTypes::Rapidity>(dY,deltaYMin(),deltaYMax(),weight) ) 
  {
    parent->lastCutWeight(0.0);
    return false;
  }
  if ( !parent->isInside<CutTypes::Azimuth>(dPhi,deltaPhiMin(),deltaPhiMax(),weight) ) 
  {
    parent->lastCutWeight(0.0);
    return false;
  }
  if ( !parent->isInside<CutTypes::Rapidity>(dR,deltaRMin(),deltaRMax(),weight) ) 
  {
    parent->lastCutWeight(0.0);
    return false;
  }

  parent->lastCutWeight(weight);
  return true;

}

void MatchboxDeltaRCut::describe() const {

  CurrentGenerator::log() 
    << fullName() << "\n"
    << "matching distances between: '"
    << theFirstMatcher->name() << "' and '"
    << theSecondMatcher->name() << "':\n"
    << "DeltaRMin = " << theDeltaRMin << " \n"
    << "DeltaRMax = " << theDeltaRMax << " \n"
    << "DeltaPhiMin = " << theDeltaPhiMin << " \n"
    << "DeltaPhiMax = " << theDeltaPhiMax << " \n"
    << "DeltaYMin = " << theDeltaYMin << " \n"
    << "DeltaYMax = " << theDeltaYMax << " \n\n";

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxDeltaRCut::persistentOutput(PersistentOStream & os) const {
  os << theDeltaYMin << theDeltaYMax 
     << theDeltaPhiMin << theDeltaPhiMax 
     << theDeltaRMin << theDeltaRMax 
     << theFirstMatcher << theSecondMatcher;
}

void MatchboxDeltaRCut::persistentInput(PersistentIStream & is, int) {
  is >> theDeltaYMin >> theDeltaYMax 
     >> theDeltaPhiMin >> theDeltaPhiMax 
     >> theDeltaRMin >> theDeltaRMax 
     >> theFirstMatcher >> theSecondMatcher;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxDeltaRCut,TwoCutBase>
  describeHerwigMatchboxDeltaRCut("Herwig::MatchboxDeltaRCut", "HwMatchboxCuts.so");

void MatchboxDeltaRCut::Init() {

  static ClassDocumentation<MatchboxDeltaRCut> documentation
    ("This class implements cuts on legoplot, rapidity and azimuthal separation, "
     "i.e. on the \\f$\\Delta R\\f$-measure and on \\f$\\Delta Y\\f$ and \\f$\\Delta \\phi\\f$. "
     "By default the cuts are only applied to coloured particles, but "
     "may optionally be applied to all particle types. ");

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaRMin
    ("DeltaRMin",
     "The minimum allowed for the legoplot distance "
     "\\f$\\Delta R_{ij}=\\sqrt{\\Delta \\phi_{ij}^2+\\Delta Y_{ij}^2}\\f$ ",
     &MatchboxDeltaRCut::theDeltaRMin, 0.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaRMax
    ("DeltaRMax",
     "The maximum allowed for the legoplot distance "
     "\\f$\\Delta R_{ij}=\\sqrt{\\Delta \\phi_{ij}^2+\\Delta Y_{ij}^2}\\f$ ",
     &MatchboxDeltaRCut::theDeltaRMax, Constants::MaxRapidity, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaPhiMin
    ("DeltaPhiMin",
     "The minimum allowed for the azimuthal separation "
     "\\f$\\Delta \\phi_{ij}\\f$ ",
     &MatchboxDeltaRCut::theDeltaPhiMin, 0.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaPhiMax
    ("DeltaPhiMax",
     "The maximum allowed for the azimuthal separation "
     "\\f$\\Delta \\phi_{ij}\\f$ ",
     &MatchboxDeltaRCut::theDeltaPhiMax, 2.0*Constants::pi, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaYMin
    ("DeltaYMin",
     "The minimum allowed for the rapidity separation "
     "\\f$\\Delta Y_{ij}\\f$ ",
     &MatchboxDeltaRCut::theDeltaYMin, 0.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaYMax
    ("DeltaYMax",
     "The maximum allowed for the rapidity separation "
     "\\f$\\Delta Y_{ij}\\f$ ",
     &MatchboxDeltaRCut::theDeltaYMax, Constants::MaxRapidity, 0, 0,
     false, false, Interface::lowerlim);

  static Reference<MatchboxDeltaRCut,MatcherBase> interfaceFirstMatcher
    ("FirstMatcher",
     "Matcher for first particle of type pitype in the pair (pitype,pjtype). "
     "If non-null only particles matching this object will be affected "
     "by the cut. ",
     &MatchboxDeltaRCut::theFirstMatcher, true, false, true, true, false);
//      &MatchboxDeltaRCut::theFirstMatcher, false, false, true, false, false);

  static Reference<MatchboxDeltaRCut,MatcherBase> interfaceSecondMatcher
    ("SecondMatcher",
     "Matcher for second particle of type pjtype in the pair (pitype,pjtype). "
     "If non-null only particles matching this object will be affected "
     "by the cut. ",
     &MatchboxDeltaRCut::theSecondMatcher, true, false, true, true, false);
//      &MatchboxDeltaRCut::theSecondMatcher, false, false, true, false, false);

}

