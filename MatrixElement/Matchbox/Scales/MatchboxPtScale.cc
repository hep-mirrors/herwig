// -*- C++ -*-
//
// MatchboxPtScale.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxPtScale class.
//

#include "MatchboxPtScale.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxPtScale::MatchboxPtScale() {}

MatchboxPtScale::~MatchboxPtScale() {}

IBPtr MatchboxPtScale::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxPtScale::fullclone() const {
  return new_ptr(*this);
}

Energy2 MatchboxPtScale::renormalizationScale() const {
  tcPDVector pd (mePartonData().begin() + 2, mePartonData().end());
  vector<LorentzMomentum> p (meMomenta().begin() + 2, meMomenta().end());
  tcPDPtr t1 = mePartonData()[0];
  tcPDPtr t2 = mePartonData()[1];
  tcCutsPtr cuts = lastCutsPtr();

  theJetFinder->cluster(pd, p, cuts, t1, t2);

  bool gotone = false;

  Energy2 maxpt2 = ZERO;
  tcPDVector::const_iterator itpd = pd.begin();
  for (vector<LorentzMomentum>::const_iterator itp = p.begin() ;
       itp != p.end(); ++itp, ++itpd )
    if ( theJetFinder->unresolvedMatcher()->check(**itpd) ) {
      gotone = true;
      maxpt2 = max(maxpt2,(*itp).perp2());
    }

  if ( !gotone && lastXCombPtr()->willPassCuts() )
    throw Exception() << "MatchboxPtScale::renormalizationScale(): No jet could be found. Check your setup."
		      << Exception::runerror;

  return maxpt2;
}

Energy2 MatchboxPtScale::factorizationScale() const {
  return renormalizationScale();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxPtScale::persistentOutput(PersistentOStream & os) const {
  os << theJetFinder;
}

void MatchboxPtScale::persistentInput(PersistentIStream & is, int) {
  is >> theJetFinder;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxPtScale,MatchboxScaleChoice>
  describeHerwigMatchboxPtScale("Herwig::MatchboxPtScale", "HwMatchboxScales.so");

void MatchboxPtScale::Init() {

  static ClassDocumentation<MatchboxPtScale> documentation
    ("MatchboxPtScale implements scale choices related to transverse momenta.");

  static Reference<MatchboxPtScale,JetFinder> interfaceJetFinder
    ("JetFinder",
     "A reference to the jet finder.",
     &MatchboxPtScale::theJetFinder, false, false, true, false, false);

}

