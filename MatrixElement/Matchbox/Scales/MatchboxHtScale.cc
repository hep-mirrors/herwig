// -*- C++ -*-
//
// MatchboxHtScale.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxHtScale class.
//

#include "MatchboxHtScale.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxHtScale::MatchboxHtScale() 
  : theIncludeMT(false), theHTFactor(1.0),
    theMTFactor(1.0) {}

MatchboxHtScale::~MatchboxHtScale() {}

IBPtr MatchboxHtScale::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxHtScale::fullclone() const {
  return new_ptr(*this);
}

Energy2 MatchboxHtScale::renormalizationScale() const {
  tcPDVector pd (mePartonData().begin() + 2, mePartonData().end());
  vector<LorentzMomentum> p (meMomenta().begin() + 2, meMomenta().end());
  tcPDPtr t1 = mePartonData()[0];
  tcPDPtr t2 = mePartonData()[1];
  tcCutsPtr cuts = lastCutsPtr();

  theJetFinder->cluster(pd, p, cuts, t1, t2);

  initWeightFactors(pd,p,theJetFinder);

  // momentum of the non-jet system
  LorentzMomentum nonJetMomentum(ZERO,ZERO,ZERO,ZERO);

  // (weighted) pt of the jet systems
  Energy ptJetSum = ZERO;

  bool gotone = false;

  tcPDVector::const_iterator pdata = pd.begin();
  vector<LorentzMomentum>::const_iterator mom = p.begin();
  for ( ; mom != p.end(); ++pdata, ++mom ) {
    if ( theJetFinder->unresolvedMatcher()->check(**pdata) ) {
      gotone = true;
      ptJetSum += jetPtWeight(*mom)*mom->perp();
    } else if ( theIncludeMT ) {
      nonJetMomentum += *mom;
    }
  }

  if ( !gotone && lastXCombPtr()->willPassCuts() )
    throw Exception() << "MatchboxHtScale::renormalizationScale(): "
		      << "No jets could be found. Check your setup."
		      << Exception::runerror;

  Energy mtNonJetSum = 
    sqrt(nonJetMomentum.perp2() + nonJetMomentum.m2());

  mtNonJetSum *= theMTFactor;
  ptJetSum *= theHTFactor;

  return sqr(ptJetSum + mtNonJetSum);

}

Energy2 MatchboxHtScale::factorizationScale() const {
  return renormalizationScale();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxHtScale::persistentOutput(PersistentOStream & os) const {
  os << theJetFinder << theIncludeMT << theHTFactor << theMTFactor;
}

void MatchboxHtScale::persistentInput(PersistentIStream & is, int) {
  is >> theJetFinder >> theIncludeMT >> theHTFactor >> theMTFactor;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxHtScale,MatchboxScaleChoice>
  describeHerwigMatchboxHtScale("Herwig::MatchboxHtScale", "HwMatchboxScales.so");

void MatchboxHtScale::Init() {

  static ClassDocumentation<MatchboxHtScale> documentation
    ("MatchboxHtScale implements scale choices related to transverse momenta.");

  static Reference<MatchboxHtScale,JetFinder> interfaceJetFinder
    ("JetFinder",
     "A reference to the jet finder.",
     &MatchboxHtScale::theJetFinder, false, false, true, false, false);

  static Switch<MatchboxHtScale,bool> interfaceIncludeMT
    ("IncludeMT",
     "Include the transverse masses of the non-jet objects.",
     &MatchboxHtScale::theIncludeMT, false, false, false);
  static SwitchOption interfaceIncludeMTYes
    (interfaceIncludeMT,
     "Yes",
     "",
     true);
  static SwitchOption interfaceIncludeMTNo
    (interfaceIncludeMT,
     "No",
     "",
     false);

  static Parameter<MatchboxHtScale,double> interfaceHTFactor
    ("HTFactor",
     "A factor to scale the HT contribution.",
     &MatchboxHtScale::theHTFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxHtScale,double> interfaceMTFactor
    ("MTFactor",
     "A factor to scale the MT contribution.",
     &MatchboxHtScale::theMTFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

}

