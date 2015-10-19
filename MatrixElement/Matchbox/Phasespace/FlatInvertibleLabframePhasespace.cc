// -*- C++ -*-
//
// FlatInvertiblePhasespace.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FlatInvertibleLabframePhasespace class.
//

#include "FlatInvertibleLabframePhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/GSLBisection.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FlatInvertibleLabframePhasespace::FlatInvertibleLabframePhasespace()
  : theLogSHat(false) {}

FlatInvertibleLabframePhasespace::~FlatInvertibleLabframePhasespace() {}

IBPtr FlatInvertibleLabframePhasespace::clone() const {
  return new_ptr(*this);
}

IBPtr FlatInvertibleLabframePhasespace::fullclone() const {
  return new_ptr(*this);
}

double FlatInvertibleLabframePhasespace::invertTwoToNKinematics(const vector<Lorentz5Momentum>& momenta,
								double* r) const {

  double weight = 1.;

  Energy finalstatemass = 0*GeV;
  for ( vector<Lorentz5Momentum>::const_iterator p =
        momenta.begin()+2; p != momenta.end(); ++p )
    finalstatemass += p->mass();

  Lorentz5Momentum pinitial = momenta[0]+momenta[1];
  Energy2 sh = pinitial.m2();
  double tau = sh/lastS();
  Energy2 shmax = lastCuts().sHatMax();
  Energy2 shmin = max(lastCuts().sHatMin(),sqr(finalstatemass));
  if (theLogSHat) {
    r[0] = log(sh/shmin)/log(shmax/shmin);
    weight *= tau*log(shmax/shmin);
  } else {
    r[0] = (sh-shmin)/(shmax-shmin);
    weight *= (shmax-shmin)/lastS();
  }
  double ltau = log(tau);
  r[1] = 0.5 - pinitial.rapidity()/ltau;
  weight *= -ltau;

  vector<Lorentz5Momentum> Pcms = momenta;
  Boost toCMS = pinitial.findBoostToCM();
  for ( vector<Lorentz5Momentum>::iterator pit =
        Pcms.begin(); pit != Pcms.end(); ++pit )
    pit->boost(toCMS);
  
  weight *= FlatInvertiblePhasespace::invertTwoToNKinematics(Pcms, r+2);

  return weight;

}


double FlatInvertibleLabframePhasespace::generateTwoToNKinematics(const double* r,
							          vector<Lorentz5Momentum>& momenta) {

  double weight = 1.;

  Energy finalstatemass = 0*GeV;
  for ( vector<Lorentz5Momentum>::const_iterator p =
        momenta.begin()+2; p != momenta.end(); ++p )
    finalstatemass += p->mass();

  Energy beamenergy = sqrt(lastS())/2.;
  Energy2 shmax = lastCuts().sHatMax();
  Energy2 shmin = max(lastCuts().sHatMin(),sqr(finalstatemass));
  Energy2 sh;
  double tau; 
  if (theLogSHat) {
    sh = shmin*pow(shmax/shmin, r[0]);
    tau = sh/lastS(); 
    weight *= tau*log(shmax/shmin);
  } else {
    sh = r[0]*(shmax-shmin)+shmin;
    tau = sh/lastS(); 
    weight *= (shmax-shmin)/lastS();
  }
  double ltau = log(tau);
  double y = ltau*(0.5 - r[1]);
  weight *= -ltau;

  double x1 = sqrt(tau)*exp(y);
  double x2 = sqrt(tau)*exp(-y);
  momenta[0] = Lorentz5Momentum(0*GeV,0*GeV,+x1*beamenergy,x1*beamenergy);
  momenta[1] = Lorentz5Momentum(0*GeV,0*GeV,-x2*beamenergy,x2*beamenergy);
  lastXCombPtr()->lastX1X2(make_pair(x1,x2));
  lastXCombPtr()->lastSHat(sh);

  weight *= FlatInvertiblePhasespace::generateTwoToNKinematics(r+2, momenta);

  // find boost to the relevant partonic frame note final state kinematics are
  // always generated in the CMS for this phase space algorithm
  Boost boostinitial = (momenta[0]+momenta[1]).findBoostToCM();
  for ( vector<Lorentz5Momentum>::iterator pit =
          momenta.begin()+2; pit != momenta.end(); ++pit )
    pit->boost(-boostinitial);

  fillDiagramWeights();

  return weight;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FlatInvertibleLabframePhasespace::persistentOutput(PersistentOStream & os) const {
  os << theLogSHat;
}

void FlatInvertibleLabframePhasespace::persistentInput(PersistentIStream & is, int) {
  is >> theLogSHat;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FlatInvertibleLabframePhasespace,MatchboxPhasespace>
  describeHerwigFlatInvertibleLabframePhasespace("Herwig::FlatInvertibleLabframePhasespace", "Herwig.so");

void FlatInvertibleLabframePhasespace::Init() {

  static ClassDocumentation<FlatInvertibleLabframePhasespace> documentation
    ("FlatInvertibleLabframePhasespace implements flat, invertible phase space generation in the lab frame.");

  static Switch<FlatInvertibleLabframePhasespace,bool> interfaceLogSHat
    ("LogSHat",
     "Generate a flat distribution in \\f$\\log(\\hat{s})\\f$.",
     &FlatInvertibleLabframePhasespace::theLogSHat, false, false, false);

  static SwitchOption interfaceLogSHatOn
    (interfaceLogSHat,
     "True", "Generate flat in \\f$\\log(\\hat{s})\\f$", true);

  static SwitchOption interfaceLogSHatOff
    (interfaceLogSHat,
     "False", "Generate flat in \\f$\\hat{s}\\f$", false);

}

