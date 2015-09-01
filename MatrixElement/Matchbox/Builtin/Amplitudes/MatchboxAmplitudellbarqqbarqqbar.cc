// -*- C++ -*-
//
// MatchboxAmplitudellbarqqbarqqbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudellbarqqbarqqbar class.
//

#include "MatchboxAmplitudellbarqqbarqqbar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudellbarqqbarqqbar::MatchboxAmplitudellbarqqbarqqbar() {}

MatchboxAmplitudellbarqqbarqqbar::~MatchboxAmplitudellbarqqbarqqbar() {}

IBPtr MatchboxAmplitudellbarqqbarqqbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudellbarqqbarqqbar::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudellbarqqbarqqbar::doinit() {
  MatchboxZGammaAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  nPoints(6);
}

void MatchboxAmplitudellbarqqbarqqbar::doinitrun() {
  MatchboxZGammaAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  nPoints(6);
}

bool MatchboxAmplitudellbarqqbarqqbar::canHandle(const PDVector& proc) const {
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
  quark = xproc.begin();
  quarkId = 0;
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
  antiQuark = xproc.begin();
  for ( ; antiQuark != xproc.end(); ++antiQuark )
    if ( (**antiQuark).id() == -quarkId ) {
      break;
    }
  if ( antiQuark == xproc.end() )
    return false;
  xproc.erase(antiQuark);
  return xproc.empty();
}

inline bool leftNonZero(int heli, int helj, int helk, int hell) {
  return 
    heli == 1 && helj == 1 &&
    abs(helk+hell) == 2;
}

inline bool rightNonZero(int heli, int helj, int helk, int hell) {
  return 
    heli == -1 && helj == -1 &&
    abs(helk+hell) == 2;
}

void MatchboxAmplitudellbarqqbarqqbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

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

Complex MatchboxAmplitudellbarqqbarqqbar::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {

  const LorentzVector<Complex>& leptonLeft
    = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& leptonRight
    = llbarRightCurrent(0,hel[0],1,hel[1]); 

  Complex LL2345 =
    leftNonZero(hel[2],hel[3],hel[4],hel[5]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    leptonLeft.dot(qqbarqqbarLeftCurrent(2,hel[2],3,hel[3],4,hel[4],5,hel[5])) : 0.;
  Complex LL4523 =
    leftNonZero(hel[4],hel[5],hel[2],hel[3]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    leptonLeft.dot(qqbarqqbarLeftCurrent(4,hel[4],5,hel[5],2,hel[2],3,hel[3])) : 0.;
  Complex LL2543 =
    leftNonZero(hel[2],hel[5],hel[4],hel[3]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    -leptonLeft.dot(qqbarqqbarLeftCurrent(2,hel[2],5,hel[5],4,hel[4],3,hel[3])) : 0.;
  Complex LL4325 =
    leftNonZero(hel[4],hel[3],hel[2],hel[5]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    -leptonLeft.dot(qqbarqqbarLeftCurrent(4,hel[4],3,hel[3],2,hel[2],5,hel[5])) : 0.;

  Complex LR2345 =
    rightNonZero(hel[2],hel[3],hel[4],hel[5]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    leptonLeft.dot(qqbarqqbarRightCurrent(2,hel[2],3,hel[3],4,hel[4],5,hel[5])) : 0.;
  Complex LR4523 =
    rightNonZero(hel[4],hel[5],hel[2],hel[3]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    leptonLeft.dot(qqbarqqbarRightCurrent(4,hel[4],5,hel[5],2,hel[2],3,hel[3])) : 0.;
  Complex LR2543 =
    rightNonZero(hel[2],hel[5],hel[4],hel[3]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    -leptonLeft.dot(qqbarqqbarRightCurrent(2,hel[2],5,hel[5],4,hel[4],3,hel[3])) : 0.;
  Complex LR4325 =
    rightNonZero(hel[4],hel[3],hel[2],hel[5]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    -leptonLeft.dot(qqbarqqbarRightCurrent(4,hel[4],3,hel[3],2,hel[2],5,hel[5])) : 0.;

  Complex RL2345 =
    leftNonZero(hel[2],hel[3],hel[4],hel[5]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    leptonRight.dot(qqbarqqbarLeftCurrent(2,hel[2],3,hel[3],4,hel[4],5,hel[5])) : 0.;
  Complex RL4523 =
    leftNonZero(hel[4],hel[5],hel[2],hel[3]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    leptonRight.dot(qqbarqqbarLeftCurrent(4,hel[4],5,hel[5],2,hel[2],3,hel[3])) : 0.;
  Complex RL2543 =
    leftNonZero(hel[2],hel[5],hel[4],hel[3]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    -leptonRight.dot(qqbarqqbarLeftCurrent(2,hel[2],5,hel[5],4,hel[4],3,hel[3])) : 0.;
  Complex RL4325 =
    leftNonZero(hel[4],hel[3],hel[2],hel[5]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    -leptonRight.dot(qqbarqqbarLeftCurrent(4,hel[4],3,hel[3],2,hel[2],5,hel[5])) : 0.;

  Complex RR2345 =
    rightNonZero(hel[2],hel[3],hel[4],hel[5]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    leptonRight.dot(qqbarqqbarRightCurrent(2,hel[2],3,hel[3],4,hel[4],5,hel[5])) : 0.;
  Complex RR4523 =
    rightNonZero(hel[4],hel[5],hel[2],hel[3]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    leptonRight.dot(qqbarqqbarRightCurrent(4,hel[4],5,hel[5],2,hel[2],3,hel[3])) : 0.;
  Complex RR2543 =
    rightNonZero(hel[2],hel[5],hel[4],hel[3]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    -leptonRight.dot(qqbarqqbarRightCurrent(2,hel[2],5,hel[5],4,hel[4],3,hel[3])) : 0.;
  Complex RR4325 =
    rightNonZero(hel[4],hel[3],hel[2],hel[5]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    -leptonRight.dot(qqbarqqbarRightCurrent(4,hel[4],3,hel[3],2,hel[2],5,hel[5])) : 0.;

  double bProp = (amplitudeMomentum(0)+amplitudeMomentum(1)).m2()/lastSHat();

  Complex gamma2345 =
    Complex(0.,-1.)*(-amplitudePartonData()[2]->iCharge()/3.)*
    (LL2345 + RL2345 + LR2345 + RR2345)/bProp;
  Complex gamma2543 =
    Complex(0.,-1.)*(-amplitudePartonData()[2]->iCharge()/3.)*
    (LL2543 + RL2543 + LR2543 + RR2543)/bProp;
  Complex gamma4523 =
    Complex(0.,-1.)*(-amplitudePartonData()[4]->iCharge()/3.)*
    (LL4523 + RL4523 + LR4523 + RR4523)/bProp;
  Complex gamma4325 =
    Complex(0.,-1.)*(-amplitudePartonData()[4]->iCharge()/3.)*
    (LL4325 + RL4325 + LR4325 + RR4325)/bProp;

  bool up2 = abs(amplitudePartonData()[2]->id()) % 2 == 0;
  bool up4 = abs(amplitudePartonData()[4]->id()) % 2 == 0;

  Complex Z2345 =
    Complex(0.,-1.)*
    (standardModel()->le()*(up2 ? standardModel()->lu() : standardModel()->ld())*LL2345 +
     standardModel()->re()*(up2 ? standardModel()->lu() : standardModel()->ld())*RL2345 +
     standardModel()->le()*(up2 ? standardModel()->ru() : standardModel()->rd())*LR2345 +
     standardModel()->re()*(up2 ? standardModel()->ru() : standardModel()->rd())*RR2345)/
    Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());
  Complex Z2543 =
    Complex(0.,-1.)*
    (standardModel()->le()*(up2 ? standardModel()->lu() : standardModel()->ld())*LL2543 +
     standardModel()->re()*(up2 ? standardModel()->lu() : standardModel()->ld())*RL2543 +
     standardModel()->le()*(up2 ? standardModel()->ru() : standardModel()->rd())*LR2543 +
     standardModel()->re()*(up2 ? standardModel()->ru() : standardModel()->rd())*RR2543)/
    Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());
  Complex Z4523 =
    Complex(0.,-1.)*
    (standardModel()->le()*(up4 ? standardModel()->lu() : standardModel()->ld())*LL4523 +
     standardModel()->re()*(up4 ? standardModel()->lu() : standardModel()->ld())*RL4523 +
     standardModel()->le()*(up4 ? standardModel()->ru() : standardModel()->rd())*LR4523 +
     standardModel()->re()*(up4 ? standardModel()->ru() : standardModel()->rd())*RR4523)/
    Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());
  Complex Z4325 =
    Complex(0.,-1.)*
    (standardModel()->le()*(up4 ? standardModel()->lu() : standardModel()->ld())*LL4325 +
     standardModel()->re()*(up4 ? standardModel()->lu() : standardModel()->ld())*RL4325 +
     standardModel()->le()*(up4 ? standardModel()->ru() : standardModel()->rd())*LR4325 +
     standardModel()->re()*(up4 ? standardModel()->ru() : standardModel()->rd())*RR4325)/
    Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());

  Complex sum2345 = 0.0;
  Complex sum2543 = 0.0;
  Complex sum4523 = 0.0;
  Complex sum4325 = 0.0;

  if ( includeGamma() ) {
    sum2345 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*gamma2345;
    sum2543 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*gamma2543;
    sum4523 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*gamma4523;
    sum4325 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*gamma4325;
  }

  if ( includeZ() ) {
    sum2345 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*Z2345;
    sum2543 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*Z2543;
    sum4523 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*Z4523;
    sum4325 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*Z4325;
  }

  double Nc = SM().Nc();

  Complex resLeading = 0.;
  Complex resSubLeading = 0.;

  if ( amplitudeToColourMap()[2] == 0 && amplitudeToColourMap()[3] == 1 &&
       amplitudeToColourMap()[4] == 2 && amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) { //(23)(45)
      resLeading = sum2543 + sum4325;
      resSubLeading = sum2345 + sum4523;
    } else if ( a == 1 ) {  //(25)(43)
      resLeading = sum2345 + sum4523;
      resSubLeading = sum2543 + sum4325;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 0 && amplitudeToColourMap()[3] == 3 &&
	      amplitudeToColourMap()[4] == 2 && amplitudeToColourMap()[5] == 1 ) {
    if ( a == 0 ) { // (25)(43)
      resLeading = sum2345 + sum4523;
      resSubLeading = sum2543 + sum4325;
    } else if ( a == 1 ) { // (23)(45)
      resLeading = sum2543 + sum4325;
      resSubLeading = sum2345 + sum4523;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 2 && amplitudeToColourMap()[3] == 3 &&
	      amplitudeToColourMap()[4] == 0 && amplitudeToColourMap()[5] == 1 ) {
    if ( a == 0 ) { //(23)(45)
      resLeading = sum2543 + sum4325;
      resSubLeading = sum2345 + sum4523;
    } else if ( a == 1 ) { //(25)(43)
      resLeading = sum2345 + sum4523;
      resSubLeading = sum2543 + sum4325;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 2 && amplitudeToColourMap()[3] == 1 &&
	      amplitudeToColourMap()[4] == 0 && amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) { //(25)(43)
      resLeading = sum2345 + sum4523;
      resSubLeading = sum2543 + sum4325;
    } else if ( a == 1 ) { //(23)(45)
      resLeading = sum2543 + sum4325;
      resSubLeading = sum2345 + sum4523;
    } else assert(false);
  } else assert(false);

  resSubLeading *= -1./Nc;
  largeN = resLeading/2.;

  return (resLeading + resSubLeading)/2.;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudellbarqqbarqqbar::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudellbarqqbarqqbar::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudellbarqqbarqqbar,MatchboxZGammaAmplitude>
  describeHerwigMatchboxAmplitudellbarqqbarqqbar("Herwig::MatchboxAmplitudellbarqqbarqqbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudellbarqqbarqqbar::Init() {

  static ClassDocumentation<MatchboxAmplitudellbarqqbarqqbar> documentation
    ("MatchboxAmplitudellbarqqbarqqbar");

}

