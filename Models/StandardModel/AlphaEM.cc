// -*- C++ -*-
//
// AlphaEM.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AlphaEM class.
//

#include "AlphaEM.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

IBPtr AlphaEM::clone() const {
  return new_ptr(*this);
}

IBPtr AlphaEM::fullclone() const {
  return new_ptr(*this);
}

void AlphaEM::persistentOutput(PersistentOStream & os) const {
  os << ounit(_me,GeV2) << ounit(_mmu,GeV2) 
     << ounit(_mtau,GeV2) << ounit(_mtop,GeV2);
}

void AlphaEM::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_me,GeV2) >> iunit(_mmu,GeV2) 
     >> iunit(_mtau,GeV2) >> iunit(_mtop,GeV2);
}

ClassDescription<AlphaEM> AlphaEM::initAlphaEM;
// Definition of the static class description member.

void AlphaEM::Init() {

  static ClassDocumentation<AlphaEM> documentation
    ("This class implements a running \\f$\\alpha_{\\mbox{EM}}\\f$ according "
     "to Buckhardt et al.",
     "In the running of $\\alpha_{EM}$, the parametrization of "
     "H.~Buckhardt et al. was used. See \\cite{KLEISSCERN9808v3pp129}.",
     "\\bibitem{KLEISSCERN9808v3pp129} R.~Kleiss et al, "
     "CERN 89-08, vol.~3, pp 129-131.");

}

void AlphaEM::doinit() {
  AlphaEMBase::doinit();
  _me   = sqr(getParticleData(ParticleID::eminus)->mass());
  _mmu  = sqr(getParticleData(ParticleID::muminus)->mass());
  _mtau = sqr(getParticleData(ParticleID::tauminus)->mass());
  _mtop = sqr(getParticleData(ParticleID::t)->mass());
}

double AlphaEM::realPi(double r) const {
  static double fvthr=1.666666666666667e0,rmax=1.e6;
  // use assymptotic formula
  if(abs(r)<1e-3) {
    return -fvthr-log(r);
  }
  // return zero for large values
  else if(abs(r)>rmax) {
    return 0.;
  }
  else if(4.*r>1.) {
    double beta=sqrt(4.*r-1.);
    return 1./3.
      -(1.+2.*r)*(2.-beta*acos(1.-1./(2.*r)));
  }
  else {
    double beta=sqrt(1.-4.*r);
    return 1./3.
      -(1.+2.*r)*(2.+beta*log(abs((beta-1.)/(beta+1.))));
  }
}

double AlphaEM::value(Energy2 scale, const StandardModelBase & sm) const {
  useMe();
  static double eps=1e-6;
  static double a1=0.0    ,b1=0.00835,c1=1.000;
  static double a2=0.0    ,b2=0.00238,c2=3.927;
  static double a3=0.00165,b3=0.00299,c3=1.000;
  static double a4=0.00221,b4=0.00293,c4=1.000;
  // alpha_EM at Q^2=0
  double alem=sm.alphaEM();
  double aempi=alem/(3.*Constants::pi);
  // convert scale to GeV2
  double Q2=scale/GeV2;
  double x=abs(Q2);
  // return q^2=0 value for small scales
  if(x<eps) return alem;
  // leptonic component
  double repigg=aempi*(realPi(_me/scale)+realPi(_mmu/scale)+realPi(_mtau/scale));
  // Hadronic component from light quarks
  if(x<9e-2)      repigg+=a1+b1*log(1.+c1*x);
  else if(x<9.)   repigg+=a2+b2*log(1.+c2*x);
  else if(x<1.e4) repigg+=a3+b3*log(1.+c3*x);
  else            repigg+=a4+b4*log(1.+c4*x);
  // Top Contribution
  repigg+=aempi*realPi(_mtop/scale);
  // reutrn the answer
  return alem/(1.-repigg);
}
