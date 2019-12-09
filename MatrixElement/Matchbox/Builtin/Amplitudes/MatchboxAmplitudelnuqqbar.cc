// -*- C++ -*-
//
// MatchboxAmplitudelnuqqbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudelnuqqbar class.
//

#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxAmplitudelnuqqbar.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;

MatchboxAmplitudelnuqqbar::MatchboxAmplitudelnuqqbar() 
        : theDiagonal(false) {}

MatchboxAmplitudelnuqqbar::~MatchboxAmplitudelnuqqbar() {}

void MatchboxAmplitudelnuqqbar::doinit() {
  MatchboxAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  theCKM = standardCKM(SM())->getUnsquaredMatrix(6);
  nPoints(4);  
}

void MatchboxAmplitudelnuqqbar::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(4);
}

IBPtr MatchboxAmplitudelnuqqbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudelnuqqbar::fullclone() const {
  return new_ptr(*this);
}

bool MatchboxAmplitudelnuqqbar::canHandle(const PDVector& proc) const {
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
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
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 && abs((*quark)->id()) % 2 == 1 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr d = *quark;
  xproc.erase(quark);
  quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 && abs((*quark)->id()) % 2 == 0 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr u = *quark;
  xproc.erase(quark);
  if ( u->iCharge() + d->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
  if ( theDiagonal && SU2Helper::family(u) != SU2Helper::family(d) ) return false;
  if ( SU2Helper::family(e) != SU2Helper::family(n) ) return false;
  return xproc.empty();
}

void MatchboxAmplitudelnuqqbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  amplitudeScale(sqrt(lastSHat()));
  setupLeptons(0,amplitudeMomentum(0),1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  MatchboxAmplitude::prepareAmplitudes(me);
}

Complex MatchboxAmplitudelnuqqbar::evaluate(size_t, const vector<int>& hel, Complex& largeN) {
  if ( abs(hel[2]+hel[3]) != 2 ) {
    largeN = 0.;
    return 0.;
  }
  Complex ckmelement = 1.;
  if ( !theDiagonal ) {
    bool wPlus = ( abs(amplitudePartonData()[0]->id()) % 2 == 0 ) ?
      amplitudePartonData()[1]->id() < 0:
      amplitudePartonData()[0]->id() < 0;
    pair<int,int> tmp(
      SU2Helper::family(amplitudePartonData()[2])-1,
      SU2Helper::family(amplitudePartonData()[3])-1);
    if ( amplitudePartonData()[3]->id() < 0 ) swap(tmp.first,tmp.second);
    ckmelement = theCKM[tmp.first][tmp.second];
    if ( !wPlus ) ckmelement = conj(ckmelement);
  }
  Complex wPropergator =
          1./Complex(((amplitudeMomentum(0)+amplitudeMomentum(1)).m2()-sqr(MW))/lastSHat(),MW*GW/lastSHat());
  Complex wVertices = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement;
  const LorentzVector<Complex>& leptonCurrent = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& quarkCurrent = qqbarLeftCurrent(2,hel[2],3,hel[3]); 
  Complex current = hel[2] == 1 ? Complex(0.,-1)*leptonCurrent.dot(quarkCurrent): 0.;
  Complex res = current*wVertices*wPropergator;
  largeN = res;
  return res;
}

Complex MatchboxAmplitudelnuqqbar::evaluateOneLoop(size_t, const vector<int>& hel) {
  if ( abs(hel[2]+hel[3]) != 2 ) return 0.;
  Complex ckmelement = 1.;
  if ( !theDiagonal ) {
    bool wPlus = ( abs(amplitudePartonData()[0]->id()) % 2 == 0 ) ?
      amplitudePartonData()[1]->id() < 0:
      amplitudePartonData()[0]->id() < 0;
    pair<int,int> tmp(
      SU2Helper::family(amplitudePartonData()[2])-1,
      SU2Helper::family(amplitudePartonData()[3])-1);
    if ( amplitudePartonData()[3]->id() < 0 ) swap(tmp.first,tmp.second);
    ckmelement = theCKM[tmp.first][tmp.second];
    if ( !wPlus ) ckmelement = conj(ckmelement);
  }
  Complex wPropergator =
          1./Complex(((amplitudeMomentum(0)+amplitudeMomentum(1)).m2()-sqr(MW))/lastSHat(),MW*GW/lastSHat());
  Complex wVertices = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement;
  const LorentzVector<Complex>& leptonCurrent = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& quarkCurrent = qqbarLeftOneLoopCurrent(2,hel[2],3,hel[3]); 
  Complex current = hel[2] == 1 ? Complex(0.,-1)*leptonCurrent.dot(quarkCurrent): 0.;
  Complex res = (SM().alphaS()/(2.*Constants::pi))*current*wVertices*wPropergator;
  return res;
}

void MatchboxAmplitudelnuqqbar::persistentOutput(PersistentOStream & os) const {
  os <<  theDiagonal << theCKM ;
}

void MatchboxAmplitudelnuqqbar::persistentInput(PersistentIStream & is, int) {
  is >>  theDiagonal >> theCKM ;
}

DescribeClass<MatchboxAmplitudelnuqqbar,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudelnuqqbar("Herwig::MatchboxAmplitudelnuqqbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudelnuqqbar::Init() {
  static ClassDocumentation<MatchboxAmplitudelnuqqbar> documentation
    ("MatchboxAmplitudelnuqqbar");
  static Switch<MatchboxAmplitudelnuqqbar,bool> interfaceDiagonal
    ("Diagonal",
     "Use a diagonal CKM matrix (ignoring the CKM object of the StandardModel).",
     &MatchboxAmplitudelnuqqbar::theDiagonal, false, false, false);
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
