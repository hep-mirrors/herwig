// -*- C++ -*-
//
// MatchboxAmplitudeggttbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudeggttbar class.
//

#include "MatchboxAmplitudeggttbar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudeggttbar::MatchboxAmplitudeggttbar() {}

MatchboxAmplitudeggttbar::~MatchboxAmplitudeggttbar() {}

IBPtr MatchboxAmplitudeggttbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudeggttbar::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudeggttbar::setupParams(map<string, double> &MGParams){
  // create parameter map for adapted Madgraph amplitude to use
  MGParams["aS"]     = SM().alphaS();
  MGParams["MZ"]     = getParticleData(ParticleID::Z0)    -> hardProcessMass()/GeV;
  MGParams["MW"]     = getParticleData(ParticleID::Wplus) -> hardProcessMass()/GeV;
  MGParams["MH"]     = getParticleData(ParticleID::h0)    -> hardProcessMass()/GeV;
  MGParams["WW"]     = getParticleData(ParticleID::Wplus) ->hardProcessWidth()/GeV;
  MGParams["WZ"]     = getParticleData(ParticleID::Z0)    ->hardProcessWidth()/GeV;
  MGParams["WH"]     = getParticleData(ParticleID::h0)    ->hardProcessWidth()/GeV;
  MGParams["MT"]     = getParticleData(ParticleID::t)     -> hardProcessMass()/GeV;
  MGParams["MB"]     = getParticleData(ParticleID::b)     -> hardProcessMass()/GeV;
  MGParams["WT"]     = getParticleData(ParticleID::t)     ->hardProcessWidth()/GeV;
  MGParams["WB"]     = getParticleData(ParticleID::b)     ->hardProcessWidth()/GeV;
  MGParams["MTA"]    = getParticleData(ParticleID::tauminus)-> hardProcessMass()/GeV;
  MGParams["GF"]     = SM().fermiConstant()*GeV2;
  MGParams["aEWM1"]  = 1./SM().alphaEMMZ();
  return;
}


void MatchboxAmplitudeggttbar::doinit() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinit();
  nPoints(4);
}

void MatchboxAmplitudeggttbar::doinitrun() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinitrun();
  nPoints(4);
}

bool MatchboxAmplitudeggttbar::canHandle(const PDVector& proc) const {
  // check process is gg > ttbar
  if ( proc.size() != 4 )
    return false;
  PDVector xproc = proc;
  PDVector::iterator top = xproc.begin();
  long topId = 0;
  for ( ; top != xproc.end(); ++top )
    if ( (**top).id() == 6 ) {
      break;
    }
  if ( top == xproc.end() )
    return false;
  topId = (**top).id();
  xproc.erase(top);
  PDVector::iterator antiTop = xproc.begin();
  for ( ; antiTop != xproc.end(); ++antiTop )
    if ( (**antiTop).id() == -topId ) {
      break;
    }
  if ( antiTop == xproc.end() )
    return false;
  xproc.erase(antiTop);
  return (xproc[0]->id()==21 && xproc[1]->id()==21);
}


void MatchboxAmplitudeggttbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  
  amplitudeScale(sqrt(lastSHat()));
  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));

  MatchboxAmplitude::prepareAmplitudes(me);
}

Complex MatchboxAmplitudeggttbar::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {

  // set up momenta to pass into Madgraph amplitude
  vector<double*> p;
  p.push_back(mom0);  p.push_back(mom1);  p.push_back(mom2);  p.push_back(mom3);
  for (int ip=0; ip<4; ++ip){
    p[ip][0] = abs(momentum(ip).e())<1.e-13 ? 0.0:double(momentum(ip).e()*amplitudeScale()/GeV);
    p[ip][1] = abs(momentum(ip).x())<1.e-13 ? 0.0:double(momentum(ip).x()*amplitudeScale()/GeV);
    p[ip][2] = abs(momentum(ip).y())<1.e-13 ? 0.0:double(momentum(ip).y()*amplitudeScale()/GeV);
    p[ip][3] = abs(momentum(ip).z())<1.e-13 ? 0.0:double(momentum(ip).z()*amplitudeScale()/GeV);
  } 

  // calculate amplitudes
  vector<complex<double> > amplitudes;
  MG_gg2ttx process;
  process.initProc(MGParams_);
  process.setMomenta(p);
  process.sigmaKin(amplitudes, hel);

  complex<double> matrixElement;
  complex<double> i = complex<double> (0, 1);

  // calculate colour flows
  if      (a==0) matrixElement = 0.;
  else if (a==1) matrixElement = i*amplitudes[0] - amplitudes[1];
  else if (a==2) matrixElement =-i*amplitudes[0] - amplitudes[2];
  else assert(false); 

  largeN = matrixElement;
  return matrixElement;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudeggttbar::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudeggttbar::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudeggttbar,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudeggttbar("Herwig::MatchboxAmplitudeggttbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudeggttbar::Init() {

  static ClassDocumentation<MatchboxAmplitudeggttbar> documentation
    ("MatchboxAmplitudeggttbar");

}

