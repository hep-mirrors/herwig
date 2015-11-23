// -*- C++ -*-
//
// MissingPtCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MissingPtCut class.
//

#include "MissingPtCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MissingPtCut::MissingPtCut() 
  : thePtMissMin(0.*GeV), thePtMissMax(Constants::MaxEnergy) {}

MissingPtCut::~MissingPtCut() {}

IBPtr MissingPtCut::clone() const {
  return new_ptr(*this);
}

IBPtr MissingPtCut::fullclone() const {
  return new_ptr(*this);
}

bool MissingPtCut::passCuts(tcCutsPtr parent, const tcPDVector & ptype, 
                            const vector<LorentzMomentum> & p) const {

  if ( thePtMissMin == ZERO && thePtMissMax == Constants::MaxEnergy )
    return true;

  // Energy ptMissSum = 0.0*GeV;
  LorentzMomentum momentumMissSum;
  bool nonu = true;

  for ( int i = 0, N = ptype.size(); i < N; ++i ) {

    if ( invisibleParticles().size() == 0 ) {
      if ( matcher()->check(*ptype[i]) ) {
        // ptMissSum = ptMissSum + p[i].perp();
        momentumMissSum = momentumMissSum + p[i];
        nonu = false;
      }
    }
    else if ( invisibleParticles().size() != 0 ) {
      for ( vector<int>::const_iterator iID = invisibleParticles().begin(); iID != invisibleParticles().end(); ++iID ) {
        int iInt = *iID;
        if ( abs(ptype[i]->id())==iInt ) {
          // ptMissSum = ptMissSum + p[i].perp();
          momentumMissSum = momentumMissSum + p[i];
          nonu = false;
        }
      }
    }

  }

  if ( nonu ) return true;

  Energy ptMiss = momentumMissSum.perp();  

  double weight = 1.0;

  // if ( !parent->isInside<CutTypes::Momentum>(ptMissSum,ptMissMin(),ptMissMax(),weight) ) {
  if ( !parent->isInside<CutTypes::Momentum>(ptMiss,ptMissMin(),ptMissMax(),weight) ) {
    parent->lastCutWeight(0.0);
    return false;
  }

  parent->lastCutWeight(weight);
  return true;

}

string MissingPtCut::doInvisibleParticles(string in) {
  istringstream ins(in);
  int first;
  ins >> first;
  theInvisibleParticles.push_back(first);
  return "";
}

void MissingPtCut::describe() const {

  CurrentGenerator::log()
    << "MissingPtCut '" << name() << "' matching "
    << "'" << matcher()->name() << "'";
  CurrentGenerator::log() << " within:\n";

  CurrentGenerator::log() 
    << "ptMiss  = " << ptMissMin()/GeV << " .. " << ptMissMax()/GeV << " GeV\n";

}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MissingPtCut::persistentOutput(PersistentOStream & os) const {
//   os << ounit(thePtMissMin,GeV) << ounit(thePtMissMax,GeV);
  os << ounit(thePtMissMin,GeV) << ounit(thePtMissMax,GeV) 
     << theInvisibleParticles << theMatcher;
}

void MissingPtCut::persistentInput(PersistentIStream & is, int) {
//   is >> iunit(thePtMissMin,GeV) >> iunit(thePtMissMax,GeV);
  is >> iunit(thePtMissMin,GeV) >> iunit(thePtMissMax,GeV)
     >> theInvisibleParticles >> theMatcher;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MissingPtCut,MultiCutBase>
  describeHerwigMissingPtCut("Herwig::MissingPtCut", "HwMatchboxCuts.so");

void MissingPtCut::Init() {

  static ClassDocumentation<MissingPtCut> documentation
//    ("MissingPtCut implements a cut on the missing transverse momentum "
//     "of a set of outgoing particles, i.e. for now the total transverse momentum "
//     "of all outgoing neutrinos in an event.");
    ("MissingPtCut implements a cut on the transverse momentum of the four-momentum "
     "sum of a set of outgoing particles that cannot be detected. By default the three "
     "standard model neutrinos are considered. If at least one undetectable particle "
     "is specified through the InvisibleParticles interface, the default choice is "
     "nullified.");

  static Command<MissingPtCut> interfaceInvisibleParticles
    ("InvisibleParticles",
     "Insert the PDG code of a particle that cannot be detected. If no particle " 
     "is inserted at all, the three standard model neutrinos are considered by "
     "default. If at least one particle is inserted, the default choice is nullified.",
     &MissingPtCut::doInvisibleParticles, false);

  static Parameter<MissingPtCut,Energy> interfacePtMissMin
    ("PtMissMin",
     "The minimum missing pt required.",
     &MissingPtCut::thePtMissMin, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<MissingPtCut,Energy> interfacePtMissMax
    ("PtMissMax",
     "The maximum missing pt allowed.",
     &MissingPtCut::thePtMissMax, GeV, Constants::MaxEnergy, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Reference<MissingPtCut,MatcherBase> interfaceMatcher
    ("Matcher",
     "A matcher for particles to cut on.",
     &MissingPtCut::theMatcher, false, false, true, false, false);

}

