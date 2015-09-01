// -*- C++ -*-
//
// MatchboxAmplitudelnuqqbargg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudelnuqqbargg class.
//

#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxAmplitudelnuqqbargg.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;

MatchboxAmplitudelnuqqbargg::MatchboxAmplitudelnuqqbargg() 
        : theDiagonal(false) {}

MatchboxAmplitudelnuqqbargg::~MatchboxAmplitudelnuqqbargg() {}

void MatchboxAmplitudelnuqqbargg::doinit() {
  MatchboxAmplitude::doinit();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  theCKM = standardCKM(SM())->getUnsquaredMatrix(6);
  nPoints(6);
}

void MatchboxAmplitudelnuqqbargg::doinitrun() {
  MatchboxAmplitude::doinitrun();
  nPoints(6);
}

IBPtr MatchboxAmplitudelnuqqbargg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudelnuqqbargg::fullclone() const {
  return new_ptr(*this);
}

bool MatchboxAmplitudelnuqqbargg::canHandle(const PDVector& proc) const {
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
  PDVector::iterator gluon = xproc.begin();
  for ( ; gluon != xproc.end(); ++gluon ) 
    if ( (*gluon)->id() == 21 ) break; 
  if ( gluon == xproc.end() ) return false;
  xproc.erase(gluon);
  gluon = xproc.begin();
  for ( ; gluon != xproc.end(); ++gluon ) 
    if ( (*gluon)->id() == 21 ) break; 
  if ( gluon == xproc.end() ) return false;
  xproc.erase(gluon);
  if ( SU2Helper::family(e) != SU2Helper::family(n) ) return false;
  if ( u->iCharge() + d->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
  if ( theDiagonal && SU2Helper::family(u) != SU2Helper::family(d) ) return false;
  return xproc.empty();
}

void MatchboxAmplitudelnuqqbargg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
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

Complex MatchboxAmplitudelnuqqbargg::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {
  if ( abs(hel[2]+hel[3]) != 2 ) {
    largeN = 0.;
    return 0.;
  }
  assert ( amplitudeToColourMap()[2] == 0 && amplitudeToColourMap()[3] == 1 );
  int g1,hg1,g2,hg2;
  if ( amplitudeToColourMap()[4] == 2 && amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) {
      g1 = 4; hg1 = hel[4];
      g2 = 5; hg2 = hel[5];
    } else if ( a == 1 ) {
      g1 = 5; hg1 = hel[5];
      g2 = 4; hg2 = hel[4];
    } else assert ( false );
  } else if ( amplitudeToColourMap()[4] == 3 && amplitudeToColourMap()[5] == 2 ) {
    if ( a == 0 ) {
      g1 = 5; hg1 = hel[5];
      g2 = 4; hg2 = hel[4];
    } else if ( a == 1 ) {
      g1 = 4; hg1 = hel[4];
      g2 = 5; hg2 = hel[5];
    } else assert ( false );
  } else assert ( false );
  Complex ckmelement = 1.;
  if ( !theDiagonal ) {
    bool wPlus = ( abs(amplitudePartonData()[0]->id()) % 2 == 0 ) ?
      amplitudePartonData()[1]->id() < 0:
      amplitudePartonData()[0]->id() < 0;
    pair<int,int> tmp(
      SU2Helper::family(amplitudePartonData()[2]),
      SU2Helper::family(amplitudePartonData()[3]));
    if ( amplitudePartonData()[3]->id() < 0 ) swap(tmp.first,tmp.second);
    ckmelement = theCKM[tmp.first][tmp.second];
    if ( !wPlus ) ckmelement = conj(ckmelement);
  }
  Complex wPropergator =
          1./Complex(((amplitudeMomentum(0)+amplitudeMomentum(1)).m2()-sqr(MW))/lastSHat(),MW*GW/lastSHat());
  Complex wVertices = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement;
  Complex sVertices =
          4.*Constants::pi*SM().alphaS();
  const LorentzVector<Complex>& leptonCurrent = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& quarkCurrent = qqbarggLeftCurrent(2,hel[2],3,hel[3],g1,hg1,g2,hg2); 
  Complex current = hel[2] == 1 ? Complex(0.,-1)*leptonCurrent.dot(quarkCurrent): 0.;
  Complex res = current*wVertices*wPropergator*sVertices;
  largeN = res;
  return res;  
}

void MatchboxAmplitudelnuqqbargg::persistentOutput(PersistentOStream & os) const {
  os << theDiagonal << theCKM << ounit(MW,GeV) << ounit(GW,GeV) << CA << CF;
}

void MatchboxAmplitudelnuqqbargg::persistentInput(PersistentIStream & is, int) {
  is >> theDiagonal >> theCKM >> iunit(MW,GeV) >> iunit(GW,GeV) >> CA >> CF;
}

DescribeClass<MatchboxAmplitudelnuqqbargg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudelnuqqbargg("Herwig::MatchboxAmplitudelnuqqbargg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudelnuqqbargg::Init() {
  static ClassDocumentation<MatchboxAmplitudelnuqqbargg> documentation
    ("MatchboxAmplitudelnuqqbargg");
  static Switch<MatchboxAmplitudelnuqqbargg,bool> interfaceDiagonal
    ("Diagonal",
     "Use a diagonal CKM matrix (ignoring the CKM object of the StandardModel).",
     &MatchboxAmplitudelnuqqbargg::theDiagonal, false, false, false);
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

