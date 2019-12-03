// -*- C++ -*-
//
// PairRapidityCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PairRapidityCut class.
//

#include "PairRapidityCut.h"
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

#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/PDT/StandardMatchers.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void PairRapidityCut::describe() const {
  CurrentGenerator::log() 
    << fullName() << "\n"
    << "matching between: '"
    << theFirstMatcher->name() << "' and '"
//     << theSecondMatcher->name() << "':\n"
    << theSecondMatcher->name() << "':\n";

//     << "y = " << theMinRapidity << " .. " << theMaxRapidity << "\n"
    for ( vector<pair<double,double> >::const_iterator r = yRanges().begin();
	  r != yRanges().end(); ++r ) {
      CurrentGenerator::log() << "y = " << r->first << " .. " << r->second << "\n";
    }

//     << "same flavour only = " << (theSameFlavourOnly?"Yes":"No") << " \n"
  CurrentGenerator::log()
    << "same flavour only = " << (theSameFlavourOnly?"Yes":"No") << " \n"
    << "opposite sign only = " << (theOppositeSignOnly?"Yes":"No") << " \n\n";
}

PairRapidityCut::PairRapidityCut() 
  : thePseudo(false), theSameFlavourOnly(false), theOppositeSignOnly(false) {}

PairRapidityCut::~PairRapidityCut() {}

IBPtr PairRapidityCut::clone() const {
  return new_ptr(*this);
}

IBPtr PairRapidityCut::fullclone() const {
  return new_ptr(*this);
}

bool PairRapidityCut::passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
		                 LorentzMomentum pi, LorentzMomentum pj,
		                 bool inci, bool incj) const {

  bool match = false;
  if ( theFirstMatcher->check(*pitype) && theSecondMatcher->check(*pjtype) ) match = true;
  if ( theFirstMatcher->check(*pjtype) && theSecondMatcher->check(*pitype) ) match = true;
//   if ( !match ||
//        ( theMinRapidity == -Constants::MaxRapidity && theMaxRapidity == Constants::MaxRapidity ) ) return true;
  if ( !match || theYRanges.empty() ) return true;
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

//   double minv = (pi+pj).rapidity();
// 
//   if ( !parent->isInside<CutTypes::Rapidity>(minv,minRapidity(),maxRapidity(),weight) ) 
//   {
//     parent->lastCutWeight(0.0);
//     return false;
//   }

//// From IdentifiedParticleCut.cc
//   double y = p.rapidity() + parent->currentYHat();
//   for ( vector<pair<double,double> >::const_iterator dy = yRanges().begin();
// 	 dy != yRanges().end(); ++dy ) {
//     if ( !parent->isInside<CutTypes::Rapidity>(y,dy->first,dy->second,weight) ) {
//       parent->lastCutWeight(0.0);
//       return false;
//     }
//   }

  double y = (pi+pj).rapidity() + parent->currentYHat();

  // Actually, why not 
  // double y = (pi+pj).rapidity() + parent->currentYHat() + parent->Y(); 
  // as in ThePEG /Cuts/Cuts.cc, /Cuts/OneCutBase.cc, etc. ???

  for ( vector<pair<double,double> >::const_iterator dy = yRanges().begin();
	dy != yRanges().end(); ++dy ) {
    if (!thePseudo) {
      if ( !parent->isInside<CutTypes::Rapidity>(y,dy->first,dy->second,weight) ) {
        parent->lastCutWeight(0.0);
        return false;
      }
    } else if (thePseudo) {
      //// From ThePEG/Cuts/OneCutBase or ThePEG/Cuts/SimpleKTCut
      // if ( p.mt()*sinh(y) <= p.perp()*sinh(theMinEta) ) return false;
      // if ( p.mt()*sinh(y) >= p.perp()*sinh(theMaxEta) ) return false;
      if ( !parent->isInside<CutTypes::Rapidity>( (pi+pj).mt()*sinh(y)/GeV ,
                                                  (pi+pj).perp()*sinh(dy->first)/GeV ,
                                                  (pi+pj).perp()*sinh(dy->second)/GeV ,
                                                  weight) ) {
        parent->lastCutWeight(0.0);
        return false;
      }
    }
  }

  parent->lastCutWeight(weight);
  return true;

}

int PairRapidityCut::family(long id) const {

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

string PairRapidityCut::doYRange(string in) {
  istringstream ins(in);
  double first, second;
  ins >> first >> second;
  if ( first > second )
    swap(first,second);
  theYRanges.push_back(make_pair(first,second));
  return "";
}

void PairRapidityCut::persistentOutput(PersistentOStream & os) const {
//   os << theMinRapidity << theMaxRapidity 
  os << theYRanges << thePseudo
     << theSameFlavourOnly << theOppositeSignOnly
     << theFirstMatcher << theSecondMatcher;
}

void PairRapidityCut::persistentInput(PersistentIStream & is, int) {
//   is >> theMinRapidity >> theMaxRapidity
  is >> theYRanges >> thePseudo
     >> theSameFlavourOnly >> theOppositeSignOnly
     >> theFirstMatcher >> theSecondMatcher;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PairRapidityCut,TwoCutBase>
  describeHerwigPairRapidityCut("Herwig::PairRapidityCut", "HwMatchboxCuts.so");

void PairRapidityCut::Init() {

  static ClassDocumentation<PairRapidityCut> documentation
    ("This class implements a rapidity cut on lepton pairs of "
     "final-state particles.");

//   static Parameter<PairRapidityCut,double> interfaceMinRapidity
//    ("MinRapidity",
//      "The minimal allowed rapditiy of the particle pair ",
//      &PairRapidityCut::theMinRapidity, -Constants::MaxRapidity, -Constants::MaxRapidity, Constants::MaxRapidity,
//      false, false, Interface::limited);
// 
//   static Parameter<PairRapidityCut,double> interfaceMaxRapidity
//    ("MaxRapidity",
//      "The maximal allowed rapidity of the particle pair ",
//      &PairRapidityCut::theMaxRapidity, Constants::MaxRapidity, -Constants::MaxRapidity, Constants::MaxRapidity,
//      false, false, Interface::limited);

  static Command<PairRapidityCut> interfaceYRange
    ("YRange",
     "Insert a rapidity range.",
     &PairRapidityCut::doYRange, false);

  static Switch<PairRapidityCut,bool> interfacePseudo
   ("Pseudo",
     "Use pseudo rapidity instead of rapidity ",
     &PairRapidityCut::thePseudo, false, false, false);
  static SwitchOption interfacePseudoNo
    (interfacePseudo,
     "No",
     "No",
     false);
  static SwitchOption interfacePseudoYes
    (interfacePseudo,
     "Yes",
     "Yes",
     true);

  static Switch<PairRapidityCut,bool> interfaceSameFlavourOnly
   ("SameFlavourOnly",
     "Whether cut works on fermion pairs of the same flavour only ",
     &PairRapidityCut::theSameFlavourOnly, true, false, false);
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

  static Switch<PairRapidityCut,bool> interfaceOppositeSignOnly
   ("OppositeSignOnly",
     "Whether cut works on fermion pairs of opposite sign only ",
     &PairRapidityCut::theOppositeSignOnly, true, false, false);
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

  static Reference<PairRapidityCut,MatcherBase> interfaceFirstMatcher
    ("FirstMatcher",
     "Matcher for first particle of type pitype in the pair (pitype,pjtype). "
     "Only particles matching this object will be affected by the cut. ",
     &PairRapidityCut::theFirstMatcher, true, false, true, true, false);

  static Reference<PairRapidityCut,MatcherBase> interfaceSecondMatcher
    ("SecondMatcher",
     "Matcher for second particle of type pjtype in the pair (pitype,pjtype). "
     "Only particles matching this object will be affected by the cut. ",
     &PairRapidityCut::theSecondMatcher, true, false, true, true, false);

}


