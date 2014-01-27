// -*- C++ -*-
//
// MatchboxHtScale.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxHtScale class.
//

#include "MatchboxHtScale.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxHtScale::MatchboxHtScale()
  : theJetsOnly(true), theDoAverage(false) {}

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

  Energy sumpt = ZERO;
  int found = 0;
  tcPDVector::const_iterator itpd = pd.begin();
  for (vector<LorentzMomentum>::const_iterator itp = p.begin() ;
       itp != p.end(); ++itp, ++itpd ) {
    if ( theJetsOnly ) {
      if ( (**itpd).coloured() ) {
	sumpt += (*itp).perp();
	found++;
      }
    }
    else {
      sumpt += (*itp).perp();
      found++;
    }
  }
  if ( theDoAverage && found != 0) return sqr(sumpt/found);

  return sqr(sumpt);
}

Energy2 MatchboxHtScale::factorizationScale() const {
  return renormalizationScale();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxHtScale::persistentOutput(PersistentOStream & os) const {
  os << theJetFinder << theJetsOnly << theDoAverage;
}

void MatchboxHtScale::persistentInput(PersistentIStream & is, int) {
  is >> theJetFinder >> theJetsOnly >> theDoAverage;
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

  static Switch<MatchboxHtScale,bool> interfaceJetsOnly
    ("JetsOnly",
     "The mode to use.",
     &MatchboxHtScale::theJetsOnly, true, false, false);
  static SwitchOption interfaceJetsOnlyTrue
    (interfaceJetsOnly,
     "True",
     "Only include jets.",
     true);
  static SwitchOption interfaceJetsOnlyFalse
    (interfaceJetsOnly,
     "False",
     "Include all outgoing particles in the matrix element.",
     false);

  static Switch<MatchboxHtScale,bool> interfaceDoAverage
    ("DoAverage",
     "The mode to use.",
     &MatchboxHtScale::theDoAverage, false, false, false);
  static SwitchOption interfaceDoAverageTrue
    (interfaceDoAverage,
     "True",
     "Average over number of jets/particles",
     true);
  static SwitchOption interfaceDoAverageFalse
    (interfaceDoAverage,
     "False",
     "Do not average.",
     false);
}

