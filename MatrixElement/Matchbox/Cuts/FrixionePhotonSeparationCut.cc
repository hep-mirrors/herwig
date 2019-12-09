// -*- C++ -*-
//
// FrixionePhotonSeparationCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FrixionePhotonSeparationCut class.
//

#include "FrixionePhotonSeparationCut.h"
#include "ThePEG/Utilities/DescribeClass.h"
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
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/PDT/StandardMatchers.h"

using namespace Herwig;

struct FrixionePartonInfo {
  double DeltaR;
  Energy pT;
  double f;
};

void FrixionePhotonSeparationCut::describe() const {
  CurrentGenerator::log() 
    << fullName()
    << " matching unresolved particles from '"
    << matcher()->name() << "':\n"
    << "DeltaZero = " << theDeltaZero << " \n"
    << "Exponent n = " << theExponentn << " \n"
    << "Efficiency = " << theEfficiency << " \n"
    << "Cut Type = " << theCutType << " \n\n";
}

IBPtr FrixionePhotonSeparationCut::clone() const {
  return new_ptr(*this);
}

IBPtr FrixionePhotonSeparationCut::fullclone() const {
  return new_ptr(*this);
}

bool FrixionePhotonSeparationCut::passCuts(tcCutsPtr parent, const tcPDVector & ptype,
					   const vector<LorentzMomentum> & p) const {
  if ( theDeltaZero <= 0 ) return true;

  double weight = 1.0;

  for ( int i = 0, N = ptype.size(); i < N; ++i ) 
    if ( ptype[i]->id() == ParticleID::gamma ) {
      vector<FrixionePartonInfo> partonvec;
      for ( int j = 0, M = ptype.size(); j < M; ++j ) {

	if ( !matcher()->check(*ptype[j]) ) continue;

        FrixionePartonInfo finfo;
        double dphi = abs(p[i].phi() - p[j].phi());
        if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
        finfo.DeltaR = sqrt(sqr(p[i].eta() - p[j].eta()) + sqr(dphi));
	if ( finfo.DeltaR < theDeltaZero ){
          finfo.pT = p[j].perp();
          finfo.f = pow((1-cos(finfo.DeltaR))/(1-cos(theDeltaZero)),theExponentn);
          partonvec.push_back(finfo);
	}

      }

      for ( unsigned int j = 0; j < partonvec.size(); ++j ) {
        Energy chidelta=ZERO;
	if (theCutType == 1) {
          for ( unsigned int k = 0; k < partonvec.size(); ++k ) 
	    if ( partonvec[k].DeltaR <= partonvec[j].DeltaR ) chidelta += partonvec[k].pT;
	}
	else if (theCutType == 2) {
          chidelta = partonvec[j].pT;
	}
        if ( !parent->isLessThan<CutTypes::Momentum>(chidelta,p[i].perp() * theEfficiency * partonvec[j].f,weight) ) {
	  parent->lastCutWeight(0.0);
	  return false;
	}
      }

    }


  parent->lastCutWeight(weight);
  return true;

}

void FrixionePhotonSeparationCut::persistentOutput(PersistentOStream & os) const {
  os << theDeltaZero << theExponentn << theEfficiency << theCutType << theMatcher;
}

void FrixionePhotonSeparationCut::persistentInput(PersistentIStream & is, int) {
  is >> theDeltaZero >> theExponentn >> theEfficiency >> theCutType >> theMatcher;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FrixionePhotonSeparationCut,MultiCutBase>
describeHerwigFrixionePhotonSeparationCut("Herwig::FrixionePhotonSeparationCut", "HwMatchboxCuts.so");

void FrixionePhotonSeparationCut::Init() {

  static ClassDocumentation<FrixionePhotonSeparationCut> documentation
    ("This class implements a separation criterium a la Frixione between "
     "final-state partons and photons.");

   static Parameter<FrixionePhotonSeparationCut,double> interfaceDeltaZero
   ("DeltaZero",
     "The maximal legoplot separation up to which partons are included in the criterium ",
     &FrixionePhotonSeparationCut::theDeltaZero, 0.7, 0.0, 10.0,
     false, false, Interface::limited);

   static Parameter<FrixionePhotonSeparationCut,double> interfaceExponentn
   ("Exponentn",
     "The exponent n of the algorithm ",
     &FrixionePhotonSeparationCut::theExponentn, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

   static Parameter<FrixionePhotonSeparationCut,double> interfaceEfficiency
   ("Efficiency",
     "The efficiency epsilon of the algorithm ",
     &FrixionePhotonSeparationCut::theEfficiency, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Switch<FrixionePhotonSeparationCut,int> interfaceCutType
    ("CutType",
     "Switch for controlling which definition of Frixione cut is used",
     &FrixionePhotonSeparationCut::theCutType, 1, false, false);
  static SwitchOption interfaceCutTypeVBFNLO
    (interfaceCutType,
     "VBFNLO",
     "Switch to Frixione cut a la VBFNLO",
     1);
  static SwitchOption interfaceCutTypeMCFM
    (interfaceCutType,
     "MCFM",
     "Switch to Frixione cut a la MCFM",
     2);

  static Reference<FrixionePhotonSeparationCut,MatcherBase> interfaceMatcher
    ("UnresolvedMatcher",
     "A matcher for particles to isolate on.",
     &FrixionePhotonSeparationCut::theMatcher, false, false, true, false, false);

  interfaceDeltaZero.setHasDefault(false);
  interfaceExponentn.setHasDefault(false);
  interfaceEfficiency.setHasDefault(false);
  interfaceCutType.setHasDefault(false);

}

