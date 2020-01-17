// -*- C++ -*-
//
// MatchboxAmplitudellbarqqbargg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudellbarqqbargg class.
//

#include "MatchboxAmplitudellbarqqbargg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudellbarqqbargg::MatchboxAmplitudellbarqqbargg() {}

MatchboxAmplitudellbarqqbargg::~MatchboxAmplitudellbarqqbargg() {}

IBPtr MatchboxAmplitudellbarqqbargg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudellbarqqbargg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudellbarqqbargg::doinit() {
  MatchboxZGammaAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  nPoints(6);
}

void MatchboxAmplitudellbarqqbargg::doinitrun() {
  MatchboxZGammaAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  nPoints(6);
}

bool MatchboxAmplitudellbarqqbargg::canHandle(const PDVector& proc) const {
  if ( proc.size() != 6 )
    return false;
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  PDVector::iterator lepton = xproc.begin();
  long leptonId = 0;
  for ( ; lepton != xproc.end(); ++lepton )
    if ( (**lepton).id() == 11 ||
	 (**lepton).id() == 13 ||
	 (**lepton).id() == 15 ) {
      break;
    }
  if ( lepton == xproc.end() )
    return false;
  leptonId = (**lepton).id();
  xproc.erase(lepton);
  PDVector::iterator antiLepton = xproc.begin();
  for ( ; antiLepton != xproc.end(); ++antiLepton )
    if ( (**antiLepton).id() == -leptonId ) {
      break;
    }
  if ( antiLepton == xproc.end() )
    return false;
  xproc.erase(antiLepton);
  PDVector::iterator quark = xproc.begin();
  long quarkId = 0;
  for ( ; quark != xproc.end(); ++quark )
    if ( abs((**quark).id()) < 6 &&
	 (**quark).id() > 0 &&
	 (**quark).hardProcessMass() == ZERO ) {
      break;
    }
  if ( quark == xproc.end() )
    return false;
  quarkId = (**quark).id();
  xproc.erase(quark);
  PDVector::iterator antiQuark = xproc.begin();
  for ( ; antiQuark != xproc.end(); ++antiQuark )
    if ( (**antiQuark).id() == -quarkId ) {
      break;
    }
  if ( antiQuark == xproc.end() )
    return false;
  xproc.erase(antiQuark);
  if ( xproc.size() != 2 )
    return false;
  return xproc[0]->id() == 21 && xproc[1]->id() == 21;
}

void MatchboxAmplitudellbarqqbargg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxZGammaAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));

  setupLeptons(0,amplitudeMomentum(0),
	       1,amplitudeMomentum(1));

  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  momentum(4,amplitudeMomentum(4));
  momentum(5,amplitudeMomentum(5));

  MatchboxZGammaAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudellbarqqbargg::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {

  if ( abs(hel[2]+hel[3]) != 2 ) {
    largeN = 0.;
    return 0.;
  }

  const LorentzVector<Complex>& leptonLeft
    = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& leptonRight
    = llbarRightCurrent(0,hel[0],1,hel[1]); 

  assert(amplitudeToColourMap()[2] == 0 &&
	 amplitudeToColourMap()[3] == 1);

  int g1,hg1,g2,hg2;
  if ( amplitudeToColourMap()[4] == 2 &&
       amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) {
      g1 = 4; hg1 = hel[4];
      g2 = 5; hg2 = hel[5];
    } else if ( a == 1 ) {
      g1 = 5; hg1 = hel[5];
      g2 = 4; hg2 = hel[4];
    } else assert(false);
  } else if ( amplitudeToColourMap()[4] == 3 &&
	      amplitudeToColourMap()[5] == 2 ) {
    if ( a == 0 ) {
      g1 = 5; hg1 = hel[5];
      g2 = 4; hg2 = hel[4];
    } else if ( a == 1 ) {
      g1 = 4; hg1 = hel[4];
      g2 = 5; hg2 = hel[5];
    } else assert(false);
  } else assert(false);

  Complex LL =
    hel[2] ==  1 ? leptonLeft.dot( qqbarggLeftCurrent (2,hel[2],3,hel[3],g1,hg1,g2,hg2))  : 0.;
  Complex RL =
    hel[2] ==  1 ? leptonRight.dot(qqbarggLeftCurrent (2,hel[2],3,hel[3],g1,hg1,g2,hg2))  : 0.;
  Complex LR =
    hel[2] == -1 ? leptonLeft.dot( qqbarggRightCurrent(2,hel[2],3,hel[3],g1,hg1,g2,hg2))  : 0.;
  Complex RR =
    hel[2] == -1 ? leptonRight.dot(qqbarggRightCurrent(2,hel[2],3,hel[3],g1,hg1,g2,hg2))  : 0.;

  double bProp = (amplitudeMomentum(0)+amplitudeMomentum(1)).m2()/lastSHat();

  Complex gamma = 0.0;
  if ( includeGamma() )
    gamma = Complex(0.,-1.)*(-amplitudePartonData()[2]->iCharge()/3.)*
      (LL + RL + LR + RR)/bProp;

  bool up = abs(amplitudePartonData()[2]->id()) % 2 == 0;
  Complex Z = 0.0;
  if ( includeZ() )
    Z = Complex(0.,-1.)*
      (standardModel()->le()*(up ? standardModel()->lu() : standardModel()->ld())*LL +
       standardModel()->re()*(up ? standardModel()->lu() : standardModel()->ld())*RL +
       standardModel()->le()*(up ? standardModel()->ru() : standardModel()->rd())*LR +
       standardModel()->re()*(up ? standardModel()->ru() : standardModel()->rd())*RR)/
      Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());

  Complex res = sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*(gamma+Z);
  largeN = res;
  return res;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudellbarqqbargg::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudellbarqqbargg::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudellbarqqbargg,MatchboxZGammaAmplitude>
  describeHerwigMatchboxAmplitudellbarqqbargg("Herwig::MatchboxAmplitudellbarqqbargg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudellbarqqbargg::Init() {

  static ClassDocumentation<MatchboxAmplitudellbarqqbargg> documentation
    ("MatchboxAmplitudellbarqqbargg");

}

