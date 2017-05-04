// -*- C++ -*-
//
// PairPtCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PairPtCut class.
//

#include "PairPtCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/PDT/StandardMatchers.h"

using namespace Herwig;

void PairPtCut::describe() const {
  CurrentGenerator::log() 
    << fullName() << "\n"
    << "matching distances between: '"
    << theFirstMatcher->name() << "' and '"
    << theSecondMatcher->name() << "':\n"
    << "pT = " << theMinPt/GeV << " .. " << theMaxPt/GeV << " GeV\n"
    << "same flavour only = " << (theSameFlavourOnly?"Yes":"No") << " \n"
    << "opposite sign only = " << (theOppositeSignOnly?"Yes":"No") << " \n\n";
}

IBPtr PairPtCut::clone() const {
  return new_ptr(*this);
}

IBPtr PairPtCut::fullclone() const {
  return new_ptr(*this);
}

bool PairPtCut::passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
		                 LorentzMomentum pi, LorentzMomentum pj,
		                 bool inci, bool incj) const {

  bool match = false;
  if ( theFirstMatcher->check(*pitype) && theSecondMatcher->check(*pjtype) ) match = true;
  if ( theFirstMatcher->check(*pjtype) && theSecondMatcher->check(*pitype) ) match = true;
  if ( !match ||
       ( theMinPt == ZERO && theMaxPt == Constants::MaxEnergy ) ) return true;
  if ( inci || incj ) return true;  

  if ( sameFlavourOnly() || oppositeSignOnly() ) {

    int fam1 = family(pitype->id());
    int fam2 = family(pjtype->id());

    if ( fam1 && fam2 ) {
      if ( sameFlavourOnly() && ( abs(fam1) != abs(fam2) ) ) return true;
      if ( oppositeSignOnly() && ( fam1*fam2 > 0 ) ) return true;
    }

  }


  double weight = 1.0;

  Energy minv = (pi+pj).perp();

  if ( !parent->isInside<CutTypes::Energy>(minv,minPt(),maxPt(),weight) ) 
  {
    parent->lastCutWeight(0.0);
    return false;
  }

  parent->lastCutWeight(weight);
  return true;

}

int PairPtCut::family(long id) const {

  int sign = (id>0)?-1:1;

  switch ( id ) {
    case ParticleID::u:
    case ParticleID::ubar:
    case ParticleID::d:
    case ParticleID::dbar:
      return 1*sign; break;
    case ParticleID::c:
    case ParticleID::cbar:
    case ParticleID::s:
    case ParticleID::sbar:
      return 2*sign; break;
    case ParticleID::t:
    case ParticleID::tbar:
    case ParticleID::b:
    case ParticleID::bbar:
      return 3*sign; break;
    case ParticleID::eminus:
    case ParticleID::eplus:
    case ParticleID::nu_e:
    case ParticleID::nu_ebar:
      return 11*sign; break;
    case ParticleID::muminus:
    case ParticleID::muplus:
    case ParticleID::nu_mu:
    case ParticleID::nu_mubar:
      return 12*sign; break;
    case ParticleID::tauminus:
    case ParticleID::tauplus:
    case ParticleID::nu_tau:
    case ParticleID::nu_taubar:
      return 13*sign; break;
  }

  return 0;

}

void PairPtCut::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMinPt,GeV) << ounit(theMaxPt,GeV) 
     << theSameFlavourOnly << theOppositeSignOnly
     << theFirstMatcher << theSecondMatcher;
}

void PairPtCut::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMinPt,GeV) >> iunit(theMaxPt,GeV) 
     >> theSameFlavourOnly >> theOppositeSignOnly
     >> theFirstMatcher >> theSecondMatcher;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PairPtCut,TwoCutBase>
  describeHerwigPairPtCut("Herwig::PairPtCut", "HwMatchboxCuts.so");

void PairPtCut::Init() {

  static ClassDocumentation<PairPtCut> documentation
    ("This class implements a transverse momentum cut on lepton pairs of "
     "final-state particles.");

  static Parameter<PairPtCut,Energy> interfaceMinPt
   ("MinPt",
     "The minimal allowed transverse momentum of the particle pair ",
     &PairPtCut::theMinPt, GeV, 0*GeV, 0*GeV, Constants::MaxEnergy,
     false, false, Interface::limited);

  static Parameter<PairPtCut,Energy> interfaceMaxPt
   ("MaxPt",
     "The maximal allowed transverse momentum of the particle pair ",
     &PairPtCut::theMaxPt, GeV, Constants::MaxEnergy, 0*GeV, Constants::MaxEnergy,
     false, false, Interface::limited);

  static Switch<PairPtCut,bool> interfaceSameFlavourOnly
   ("SameFlavourOnly",
     "Whether cut works on fermion pairs of the same flavour only ",
     &PairPtCut::theSameFlavourOnly, true, false, false);
  static SwitchOption interfaceSameFlavourOnlyYes
    (interfaceSameFlavourOnly,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceSameFlavourOnlyNo
    (interfaceSameFlavourOnly,
     "No",
     "No",
     false);

  static Switch<PairPtCut,bool> interfaceOppositeSignOnly
   ("OppositeSignOnly",
     "Whether cut works on fermion pairs of opposite sign only ",
     &PairPtCut::theOppositeSignOnly, true, false, false);
  static SwitchOption interfaceOppositeSignOnlyYes
    (interfaceOppositeSignOnly,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceOppositeSignOnlyNo
    (interfaceOppositeSignOnly,
     "No",
     "No",
     false);

  static Reference<PairPtCut,MatcherBase> interfaceFirstMatcher
    ("FirstMatcher",
     "Matcher for first particle of type pitype in the pair (pitype,pjtype). "
     "Only particles matching this object will be affected by the cut. ",
     &PairPtCut::theFirstMatcher, true, false, true, true, false);

  static Reference<PairPtCut,MatcherBase> interfaceSecondMatcher
    ("SecondMatcher",
     "Matcher for second particle of type pjtype in the pair (pitype,pjtype). "
     "Only particles matching this object will be affected by the cut. ",
     &PairPtCut::theSecondMatcher, true, false, true, true, false);

}


