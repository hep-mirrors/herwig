// -*- C++ -*-
//
// MatchboxAmplitudeqqbarttbarg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudeqqbarttbarg class.
//

#include "MatchboxAmplitudeqqbarttbarg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudeqqbarttbarg::MatchboxAmplitudeqqbarttbarg() {}

MatchboxAmplitudeqqbarttbarg::~MatchboxAmplitudeqqbarttbarg() {}

IBPtr MatchboxAmplitudeqqbarttbarg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudeqqbarttbarg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudeqqbarttbarg::setupParams(map<string, double> &MGParams){
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


void MatchboxAmplitudeqqbarttbarg::doinit() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinit();
  nPoints(5);
}

void MatchboxAmplitudeqqbarttbarg::doinitrun() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinitrun();
  nPoints(5);
}

bool MatchboxAmplitudeqqbarttbarg::canHandle(const PDVector& proc) const {
  // check process is qqbar > ttbarg or some crossing of this
  if ( proc.size() != 5 )
    return false;
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();  
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
  PDVector::iterator quark = xproc.begin();
  long quarkId = 0;
  for ( ; quark != xproc.end(); ++quark )
    if ( abs((**quark).id()) < 6 &&
	 (**quark).id() > 0  ) {
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
  return xproc[0]->id() == 21;
}


void MatchboxAmplitudeqqbarttbarg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  amplitudeScale(sqrt(lastSHat()));
  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  momentum(4,amplitudeMomentum(4));

  MatchboxAmplitude::prepareAmplitudes(me);
}

Complex MatchboxAmplitudeqqbarttbarg::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {

  // check if doing qqbar, qg or qbarg initiated process
  int crossed;
  cPDVector partons = mePartonData();
  if (abs(partons[0]->id())< 6 && 
      partons[1]->id()==-partons[0]->id()) crossed=1;
  else if (((partons[0]->id()>0 && partons[0]->id()<6) && 
  	    partons[1]->id()==ParticleID::g) ||
  	   ((partons[1]->id()>0 && partons[1]->id()<6) && 
  	    partons[0]->id()==ParticleID::g)) crossed=2;
  else if (((partons[0]->id()<0 && partons[0]->id()>-6) && 
  	    partons[1]->id()==ParticleID::g) ||
  	   ((partons[1]->id()<0 && partons[1]->id()>-6) && 
  	    partons[0]->id()==ParticleID::g)) crossed=3;
  else 
    throw Exception() << "MatchboxAmplitudeqqbarststbarg::evaluate(): Unrecognised process\n"
  		      << Exception::runerror;


  // set up momenta to pass into Madgraph amplitude
  vector<double*> p;
  p.push_back(mom0);  p.push_back(mom1);  p.push_back(mom2);  p.push_back(mom3); p.push_back(mom4);
  for (int ip=0; ip<5; ++ip){
    p[ip][0] = abs(momentum(ip).e())<1.e-13 ? 0.0:double(momentum(ip).e()*amplitudeScale()/GeV);
    p[ip][1] = abs(momentum(ip).x())<1.e-13 ? 0.0:double(momentum(ip).x()*amplitudeScale()/GeV);
    p[ip][2] = abs(momentum(ip).y())<1.e-13 ? 0.0:double(momentum(ip).y()*amplitudeScale()/GeV);
    p[ip][3] = abs(momentum(ip).z())<1.e-13 ? 0.0:double(momentum(ip).z()*amplitudeScale()/GeV);
  } 
 
  // calculate amplitudes
  vector<complex<double> > amplitudes;
  MG_qqx2ttxg process;
  double factor=lastSHat()/GeV2;
  process.initProc(MGParams_);  
  process.setMomenta(p);
  process.sigmaKin(amplitudes, hel, crossed);

  for (int iamp=0; iamp<int(amplitudes.size()); ++iamp)
    amplitudes[iamp]*=sqrt(factor);

  complex<double> matrixElement;
  complex<double> i = complex<double>(0,1);

  // calculate colour flows
  if (a==0)
    matrixElement = (1./6.)*(amplitudes[1] + amplitudes[2]);
  else if (a==1)
    matrixElement =  1./2.*( i*amplitudes[0] - amplitudes[1] - amplitudes[3]);  
  else if (a==2) 
    matrixElement = (1./6.)*( amplitudes[3] + amplitudes[4]);
  else if (a==3)
    matrixElement =  1./2.*(-i*amplitudes[0] - amplitudes[2] - amplitudes[4]);
  else assert(false);

  largeN = matrixElement;
  return matrixElement;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudeqqbarttbarg::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudeqqbarttbarg::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudeqqbarttbarg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudeqqbarttbarg("Herwig::MatchboxAmplitudeqqbarttbarg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudeqqbarttbarg::Init() {

  static ClassDocumentation<MatchboxAmplitudeqqbarttbarg> documentation
    ("MatchboxAmplitudeqqbarttbarg");

}

