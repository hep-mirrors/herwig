// -*- C++ -*-
//
// GenericOLPAmplitude.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericOLPAmplitude class.
//

#include "GenericOLPAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

GenericOLPAmplitude::GenericOLPAmplitude() {}

GenericOLPAmplitude::~GenericOLPAmplitude() {}

IBPtr GenericOLPAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr GenericOLPAmplitude::fullclone() const {
  return new_ptr(*this);
}

void GenericOLPAmplitude::signOLP(const string&, const string&) {
}

void GenericOLPAmplitude::startOLP(const string&, int& status) {
  status = 1;
}

LorentzVector<Complex> GenericOLPAmplitude::plusPolarization(const Lorentz5Momentum&,
							     const Lorentz5Momentum&,
							     int) const {

  return LorentzVector<Complex>(0.0,0.0,0.0,0.0);

}

void GenericOLPAmplitude::evalSubProcess() const {

  double units = pow(lastSHat()/GeV2,mePartonData().size()-4.);
  fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
  //double as = SM().alphaS();
  //double scale = sqrt(mu2()/GeV2);

  lastTreeME2(0.0*units);
  lastOneLoopInterference(0.0*units);
  lastOneLoopPoles(pair<double,double>(0.0*units,0.0*units));

}

void GenericOLPAmplitude::evalColourCorrelator(pair<int,int> ij) const {

  double units = pow(lastSHat()/GeV2,mePartonData().size()-4.);
  fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
  //double as = SM().alphaS();
  //double scale = sqrt(mu2()/GeV2);

  lastColourCorrelator(ij,0.0*units);

}

void GenericOLPAmplitude::evalSpinColourCorrelator(pair<int,int> ij) const {

  double units = pow(lastSHat()/GeV2,mePartonData().size()-4.);
  fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
  //double as = SM().alphaS();
  //double scale = sqrt(mu2()/GeV2);

  lastColourSpinCorrelator(ij,0.0*units);

}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void GenericOLPAmplitude::persistentOutput(PersistentOStream &) const {}

void GenericOLPAmplitude::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<GenericOLPAmplitude,MatchboxOLPME>
  describeHerwigGenericOLPAmplitude("Herwig::GenericOLPAmplitude", "Herwig.so");

void GenericOLPAmplitude::Init() {

  static ClassDocumentation<GenericOLPAmplitude> documentation
    ("GenericOLPAmplitude implements OLP interfaces.");

}

