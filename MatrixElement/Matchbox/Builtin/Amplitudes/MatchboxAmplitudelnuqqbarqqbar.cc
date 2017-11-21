// -*- C++ -*-
//
// MatchboxAmplitudelnuqqbarqqbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudelnuqqbarqqbar class.
//

#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxAmplitudelnuqqbarqqbar.h"

using namespace Herwig;

MatchboxAmplitudelnuqqbarqqbar::MatchboxAmplitudelnuqqbarqqbar() 
        : theDiagonal(false) {}

MatchboxAmplitudelnuqqbarqqbar::~MatchboxAmplitudelnuqqbarqqbar() {}

void MatchboxAmplitudelnuqqbarqqbar::doinit() {
  MatchboxAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  theCKM = standardCKM(SM())->getUnsquaredMatrix(6);
  nPoints(6);
}

void MatchboxAmplitudelnuqqbarqqbar::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(6);
}

IBPtr MatchboxAmplitudelnuqqbarqqbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudelnuqqbarqqbar::fullclone() const {
  return new_ptr(*this);
}

bool MatchboxAmplitudelnuqqbarqqbar::canHandle(const PDVector& proc) const {
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  //  Charge charge(ZERO);
  PDVector::iterator elektron = xproc.begin();
  for ( ; elektron != xproc.end(); ++elektron ) 
    if ( abs((*elektron)->id()) >= 11 && abs((*elektron)->id()) <= 16 && abs((*elektron)->id()) % 2 == 1 ) break;
  if ( elektron == xproc.end() ) return false;
  PDPtr e = *elektron;
  xproc.erase(elektron);
  PDVector::iterator neutrino = xproc.begin();
  for ( ; neutrino != xproc.end(); ++neutrino ) 
    if ( abs((*neutrino)->id()) >= 11 && abs((*neutrino)->id()) <= 16 && abs((*neutrino)->id()) % 2 == 0 ) break;
  if ( neutrino == xproc.end() ) return false;
  PDPtr n = *neutrino;
  xproc.erase(neutrino);
  PDVector::iterator quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr q1 = *quark;
  xproc.erase(quark);
  quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr q2 = *quark;
  xproc.erase(quark);
  quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr q3 = *quark;
  xproc.erase(quark);
  quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr q4 = *quark;
  xproc.erase(quark);
  if ( q1->id() == -q2->id() ) {
    if ( q3->iCharge() + q4->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q3) != SU2Helper::family(q4) ) return false;
  } else if ( q1->id() == -q3->id() ) {
    if ( q2->iCharge() + q4->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q2) != SU2Helper::family(q4) ) return false;
  } else  if ( q1->id() == -q4->id() ) {
    if ( q2->iCharge() + q3->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q2) != SU2Helper::family(q3) ) return false;
  } else  if ( q2->id() == -q3->id() ) {
    if ( q1->iCharge() + q4->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q1) != SU2Helper::family(q4) ) return false;
  } else  if ( q2->id() == -q4->id() ) {
    if ( q1->iCharge() + q3->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q1) != SU2Helper::family(q3) ) return false;
  } else  if ( q3->id() == -q4->id() ) {
    if ( q1->iCharge() + q2->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q1) != SU2Helper::family(q2) ) return false;
  } else return false;
  if ( SU2Helper::family(e) != SU2Helper::family(n) ) return false;
  return xproc.empty();
}

void MatchboxAmplitudelnuqqbarqqbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  amplitudeScale(sqrt(lastSHat()));
  setupLeptons(0,amplitudeMomentum(0),1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  momentum(4,amplitudeMomentum(4));
  momentum(5,amplitudeMomentum(5));
  MatchboxAmplitude::prepareAmplitudes(me);
}

Complex MatchboxAmplitudelnuqqbarqqbar::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {
  const LorentzVector<Complex>& leptonCurrent = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  Complex Current2345 =
    hel[2] == 1 && hel[3] == 1 && abs(hel[4]+hel[5]) == 2 &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    leptonCurrent.dot(qqbarqqbarLeftCurrent(2,hel[2],3,hel[3],4,hel[4],5,hel[5])) : 0.;
  Complex Current4523 =
    hel[4] == 1 && hel[5] == 1 && abs(hel[2]+hel[3]) == 2 &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    leptonCurrent.dot(qqbarqqbarLeftCurrent(4,hel[4],5,hel[5],2,hel[2],3,hel[3])) : 0.;
  Complex Current2543 =
    hel[2] == 1 && hel[5] == 1 && abs(hel[4]+hel[3]) == 2 &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    -leptonCurrent.dot(qqbarqqbarLeftCurrent(2,hel[2],5,hel[5],4,hel[4],3,hel[3])) : 0.;
  Complex Current4325 =
    hel[4] == 1 && hel[3] == 1 && abs(hel[2]+hel[5]) == 2 &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    -leptonCurrent.dot(qqbarqqbarLeftCurrent(4,hel[4],3,hel[3],2,hel[2],5,hel[5])) : 0.;
  Complex ckmelement23 = 1.;
  Complex ckmelement45 = 1.;
  if ( !theDiagonal ) {
    bool wPlus = ( abs(amplitudePartonData()[0]->id()) % 2 == 0 ) ?
      amplitudePartonData()[1]->id() < 0:
      amplitudePartonData()[0]->id() < 0;
    pair<int,int> tmp23(
      SU2Helper::family(amplitudePartonData()[2])-1,
      SU2Helper::family(amplitudePartonData()[3])-1);
    pair<int,int> tmp45(
      SU2Helper::family(amplitudePartonData()[4])-1,
      SU2Helper::family(amplitudePartonData()[5])-1);
    if ( amplitudePartonData()[3]->id() < 0 ) swap(tmp23.first,tmp23.second);
    if ( amplitudePartonData()[5]->id() < 0 ) swap(tmp45.first,tmp45.second);
    ckmelement23 = theCKM[tmp23.first][tmp23.second];
    ckmelement45 = theCKM[tmp45.first][tmp45.second];
    if ( !wPlus ) {
      ckmelement23 = conj(ckmelement23);
      ckmelement45 = conj(ckmelement45);
    }
  }
  Complex wPropergator =
          1./Complex(((amplitudeMomentum(0)+amplitudeMomentum(1)).m2()-sqr(MW))/lastSHat(),MW*GW/lastSHat());
  Complex wVertices23 = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement23;
  Complex wVertices45 = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement45;
  Complex sVertices =
          4.*Constants::pi*SM().alphaS();
  Complex res2345 =
    Complex(0.,-1.)*wPropergator*sVertices*Current2345*wVertices23;
  Complex res2543 =
    Complex(0.,-1.)*wPropergator*sVertices*Current2543*wVertices23;
  Complex res4523 =
    Complex(0.,-1.)*wPropergator*sVertices*Current4523*wVertices45;
  Complex res4325 =
    Complex(0.,-1.)*wPropergator*sVertices*Current4325*wVertices45;
  double Nc = SM().Nc();
  Complex resLeading = 0.;
  Complex resSubLeading = 0.;
  if ( amplitudeToColourMap()[2] == 0 && amplitudeToColourMap()[3] == 1 &&
       amplitudeToColourMap()[4] == 2 && amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) { //(23)(45)
      resLeading = res2543 + res4325;
      resSubLeading = res2345 + res4523;
    } else if ( a == 1 ) {  //(25)(43)
      resLeading = res2345 + res4523;
      resSubLeading = res2543 + res4325;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 0 && amplitudeToColourMap()[3] == 3 &&
	      amplitudeToColourMap()[4] == 2 && amplitudeToColourMap()[5] == 1 ) {
    if ( a == 0 ) { // (25)(43)
      resLeading = res2345 + res4523;
      resSubLeading = res2543 + res4325;
    } else if ( a == 1 ) { // (23)(45)
      resLeading = res2543 + res4325;
      resSubLeading = res2345 + res4523;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 2 && amplitudeToColourMap()[3] == 3 &&
	      amplitudeToColourMap()[4] == 0 && amplitudeToColourMap()[5] == 1 ) {
    if ( a == 0 ) { //(23)(45)
      resLeading = res2543 + res4325;
      resSubLeading = res2345 + res4523;
    } else if ( a == 1 ) { //(25)(43)
      resLeading = res2345 + res4523;
      resSubLeading = res2543 + res4325;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 2 && amplitudeToColourMap()[3] == 1 &&
	      amplitudeToColourMap()[4] == 0 && amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) { //(25)(43)
      resLeading = res2345 + res4523;
      resSubLeading = res2543 + res4325;
    } else if ( a == 1 ) { //(23)(45)
      resLeading = res2543 + res4325;
      resSubLeading = res2345 + res4523;
    } else assert(false);
  } else assert(false);
  resSubLeading *= -1./Nc;
  largeN = resLeading/2.;
  return (resLeading + resSubLeading)/2.;
}

void MatchboxAmplitudelnuqqbarqqbar::persistentOutput(PersistentOStream & os) const {
  os << theDiagonal << theCKM ;
}

void MatchboxAmplitudelnuqqbarqqbar::persistentInput(PersistentIStream & is, int) {
  is >> theDiagonal >> theCKM ;
}

DescribeClass<MatchboxAmplitudelnuqqbarqqbar,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudelnuqqbarqqbar("Herwig::MatchboxAmplitudelnuqqbarqqbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudelnuqqbarqqbar::Init() {
  static ClassDocumentation<MatchboxAmplitudelnuqqbarqqbar> documentation
    ("MatchboxAmplitudelnuqqbarqqbar");
  static Switch<MatchboxAmplitudelnuqqbarqqbar,bool> interfaceDiagonal
    ("Diagonal",
     "Use a diagonal CKM matrix (ignoring the CKM object of the StandardModel).",
     &MatchboxAmplitudelnuqqbarqqbar::theDiagonal, false, false, false);
  static SwitchOption interfaceDiagonalYes
    (interfaceDiagonal,
     "Yes",
     "Use a diagonal CKM matrix.",
     true);
  static SwitchOption interfaceDiagonalNo
    (interfaceDiagonal,
     "No",
     "Use the CKM object as used by the StandardModel.",
     false);
}

