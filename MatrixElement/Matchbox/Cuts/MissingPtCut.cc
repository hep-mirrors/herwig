// -*- C++ -*-
//
// MissingPtCut.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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

  Energy ptMissSum = 0.0*GeV;
  bool nonu = false;

  for ( int i = 0, N = ptype.size(); i < N; ++i ) {
    if ( abs(ptype[i]->id())==ParticleID::nu_e || abs(ptype[i]->id())==ParticleID::nu_mu || abs(ptype[i]->id())==ParticleID::nu_tau ) {
      ptMissSum = ptMissSum + p[i].perp();
      nonu = true;
    }
  }

  if ( !nonu ) return true;

  double weight = 1.0;

  if ( !parent->isInside<CutTypes::Momentum>(ptMissSum,ptMissMin(),ptMissMax(),weight) ) {
    parent->lastCutWeight(0.0);
    return false;
  }

  parent->lastCutWeight(weight);
  return true;

}

void MissingPtCut::describe() const {

  CurrentGenerator::log()
    << "MissingPtCut '" << name() << "' matching ";
  CurrentGenerator::log() << " within:\n";

  CurrentGenerator::log() 
    << "ptMiss  = " << ptMissMin()/GeV << " .. " << ptMissMax()/GeV << " GeV\n";

}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MissingPtCut::persistentOutput(PersistentOStream & os) const {
  os << ounit(thePtMissMin,GeV) << ounit(thePtMissMax,GeV);
}

void MissingPtCut::persistentInput(PersistentIStream & is, int) {
  is >> iunit(thePtMissMin,GeV) >> iunit(thePtMissMax,GeV);
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
    ("MissingPtCut implements a cut on the total missing transverse momentum, "
     "i.e. for now the total transverse momentum of all neutrinos.");

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

}

