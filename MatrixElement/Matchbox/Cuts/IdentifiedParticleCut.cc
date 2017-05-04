// -*- C++ -*-
//
// IdentifiedParticleCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IdentifiedParticleCut class.
//

#include "IdentifiedParticleCut.h"
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

IdentifiedParticleCut::IdentifiedParticleCut() 
  : thePtMin(0.*GeV), thePtMax(Constants::MaxEnergy) {}

IdentifiedParticleCut::~IdentifiedParticleCut() {}

IBPtr IdentifiedParticleCut::clone() const {
  return new_ptr(*this);
}

IBPtr IdentifiedParticleCut::fullclone() const {
  return new_ptr(*this);
}

bool IdentifiedParticleCut::passCuts(tcCutsPtr parent,
				     tcPDPtr ptype, LorentzMomentum p) const {

  if ( !matcher()->check(*ptype) ||
       ( thePtMin == ZERO && thePtMax == Constants::MaxEnergy &&
	 theYRanges.empty() ) )
    return true;

  double weight = 1.0;

  if ( !parent->isInside<CutTypes::Momentum>(p.perp(),ptMin(),ptMax(),weight) ) {
    parent->lastCutWeight(0.0);
    return false;
  }

  double y = p.rapidity() + parent->currentYHat();
  for ( vector<pair<double,double> >::const_iterator dy = yRanges().begin();
	dy != yRanges().end(); ++dy ) {
    if ( !parent->isInside<CutTypes::Rapidity>(y,dy->first,dy->second,weight) ) {
      parent->lastCutWeight(0.0);
      return false;
    }
  }

  parent->lastCutWeight(weight);
  return true;

}

void IdentifiedParticleCut::describe() const {

  CurrentGenerator::log()
    << "IdentifiedParticleCut '" << name() << "' matching "
    << "'" << matcher()->name() << "'";
  CurrentGenerator::log() << " within:\n";

  CurrentGenerator::log() 
    << "pt  = " << ptMin()/GeV << " .. " << ptMax()/GeV << " GeV\n";

  for ( vector<pair<double,double> >::const_iterator r = yRanges().begin();
	r != yRanges().end(); ++r ) {
    CurrentGenerator::log() << "y   = " << r->first << " .. " << r->second << "\n";
  }

}

string IdentifiedParticleCut::doYRange(string in) {
  istringstream ins(in);
  double first, second;
  ins >> first >> second;
  if ( first > second )
    swap(first,second);
  theYRanges.push_back(make_pair(first,second));
  return "";
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IdentifiedParticleCut::persistentOutput(PersistentOStream & os) const {
  os << ounit(thePtMin,GeV) << ounit(thePtMax,GeV)
     << theYRanges << theMatcher;
}

void IdentifiedParticleCut::persistentInput(PersistentIStream & is, int) {
  is >> iunit(thePtMin,GeV) >> iunit(thePtMax,GeV)
     >> theYRanges >> theMatcher;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IdentifiedParticleCut,OneCutBase>
  describeHerwigIdentifiedParticleCut("Herwig::IdentifiedParticleCut", "HwMatchboxCuts.so");

void IdentifiedParticleCut::Init() {

  static ClassDocumentation<IdentifiedParticleCut> documentation
    ("IdentifiedParticleCut implements cuts on single momenta.");

  static Parameter<IdentifiedParticleCut,Energy> interfacePtMin
    ("PtMin",
     "The minimum pt required.",
     &IdentifiedParticleCut::thePtMin, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<IdentifiedParticleCut,Energy> interfacePtMax
    ("PtMax",
     "The maximum pt allowed.",
     &IdentifiedParticleCut::thePtMax, GeV, Constants::MaxEnergy, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Command<IdentifiedParticleCut> interfaceYRange
    ("YRange",
     "Insert a rapidity range.",
     &IdentifiedParticleCut::doYRange, false);

  static Reference<IdentifiedParticleCut,MatcherBase> interfaceMatcher
    ("Matcher",
     "A matcher for particles to cut on.",
     &IdentifiedParticleCut::theMatcher, false, false, true, false, false);

}

