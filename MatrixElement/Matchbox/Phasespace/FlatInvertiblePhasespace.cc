// -*- C++ -*-
//
// FlatInvertiblePhasespace.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FlatInvertiblePhasespace class.
//

#include "FlatInvertiblePhasespace.h"
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

FlatInvertiblePhasespace::FlatInvertiblePhasespace() {}

FlatInvertiblePhasespace::~FlatInvertiblePhasespace() {}

IBPtr FlatInvertiblePhasespace::clone() const {
  return new_ptr(*this);
}

IBPtr FlatInvertiblePhasespace::fullclone() const {
  return new_ptr(*this);
}

double FlatInvertiblePhasespace::bisect(double v, double n, 
					double target, double maxLevel) const {

  if ( v != 0.0 && v != 1.0 ) {

    double level = 0;
    double left = 0;
    double right = 1;

    double checkV = -1.;
    double u = -1;

    while ( level < maxLevel ) {

      u = (left+right)*pow(0.5,level+1.);
      checkV = 
	pow(u,n+1.)*(n+2.-(n+1.)*u);

      if ( log10(abs(1.-checkV/v)) <= target )
	break;

      left *= 2.;
      right *= 2.;

      if ( v <= checkV ) {
	right -= 1.;
	++level;
      }

      if ( v > checkV ) {
	left += 1.;
	++level;
      }

    }

    return u;

  }

  return v;

}

static double flatWeights[7] = {

  -1.,-1.,
  0.039788735772973833942,
  0.00012598255637968550463,
  1.3296564302788840628E-7,
  7.0167897579949011130E-11,
  2.2217170114046130768E-14

};

double FlatInvertiblePhasespace::generateIntermediates(vector<Energy>& K,
						       const double* r) const {

  size_t n = K.size() + 1;
  for ( size_t i = 2; i <= n-1; ++i ) {
    double u = bisect(r[i-2],n-1-i);
    K[i-1] = sqrt(u*sqr(K[i-2]));
  } 

  return flatWeights[n];

}

double FlatInvertiblePhasespace::invertIntermediates(const vector<Energy>& K,
						     double* r) const {

  size_t n = K.size() + 1;
  for ( size_t i = 2; i <= n-1; ++i ) {
    double u = sqr(K[i-1]/K[i-2]);
    r[i-2] = (n+1-i)*pow(u,(double)(n-i)) - (n-i)*pow(u,(double)(n+1-i));
  } 

  return flatWeights[n];

}

double FlatInvertiblePhasespace::generateIntermediates(vector<Energy>& M,
						       const vector<Energy>& m,
						       const double* r) const {

  size_t n = M.size() + 1;

  vector<Energy> K = M;
  for ( size_t i = 1; i <= n; ++i )
    K[0] -= m[i-1];

  double w0 = generateIntermediates(K,r);

  M = K;
  for ( size_t i = 1; i <= n-1; ++i ) {
    for ( size_t k = i; k <= n; ++k )
      M[i-1] += m[k-1];
  }

  double weight = 8.*w0*rho(M[n-2],m[n-1],m[n-2]);

  for ( size_t i = 2; i <= n-1; ++i ) {
    weight *= 
      (rho(M[i-2],M[i-1],m[i-2])/rho(K[i-2],K[i-1],ZERO)) * (M[i-1]/K[i-1]);
  }

  weight *= pow(K[0]/M[0],2.*n-4.);

  return weight;

}

double FlatInvertiblePhasespace::invertIntermediates(const vector<Energy>& M,
						     const vector<Energy>& m,
						     double* r) const {

  size_t n = M.size() + 1;

  vector<Energy> K = M;
  for ( size_t i = 1; i <= n-1; ++i ) {
    for ( size_t k = i; k <= n; ++k )
      K[i-1] -= m[k-1];
  }

  double w0 = invertIntermediates(K,r);

  double weight = 8.*w0*rho(M[n-2],m[n-1],m[n-2]);

  for ( size_t i = 2; i <= n-1; ++i ) {
    weight *= 
      (rho(M[i-2],M[i-1],m[i-2])/rho(K[i-2],K[i-1],ZERO)) * (M[i-1]/K[i-1]);
  }

  weight *= pow(K[0]/M[0],2.*n-4.);

  return weight;

}


double FlatInvertiblePhasespace::generateKinematics(vector<Lorentz5Momentum>& P,
						    Energy Ecm,
						    const double* r) const {

  vector<Energy> m;
  for ( vector<Lorentz5Momentum>::const_iterator p =
	  P.begin() + 2; p != P.end(); ++p )
    m.push_back(p->mass());

  size_t n = P.size() - 2;
  vector<Energy> M(n-1);
  M[0] = Ecm;

  double weight = generateIntermediates(M,m,r);

  M.push_back(m.back());

  Lorentz5Momentum Q(M[0]);
  Lorentz5Momentum nextQ;

  for ( size_t i = 2; i <= n; ++i ) {

    Energy q = 4.*M[i-2]*rho(M[i-2],M[i-1],m[i-2]);

    double c = 2.*r[n-6+2*i]-1.;
    double s = sqrt(1.-sqr(c));
    double phi = 2.*Constants::pi*r[n-5+2*i];
    double cphi = cos(phi);
    double sphi = sqrt(1.-sqr(cphi));
    if ( phi > Constants::pi )
      sphi = -sphi;

    P[i].setX(q*cphi*s);
    P[i].setY(q*sphi*s);
    P[i].setZ(q*c);
    P[i].rescaleEnergy();
    P[i].boost(Q.boostVector());
    P[i].rescaleEnergy();

    nextQ = Q - P[i];
    nextQ.setMass(M[i-1]);
    nextQ.rescaleEnergy();

    Q = nextQ;

  }

  P.back() = Q;

  return weight;

}

double FlatInvertiblePhasespace::invertKinematics(const vector<Lorentz5Momentum>& P,
						  Energy Ecm,
						  double* r) const {

  vector<Energy> m;
  for ( vector<Lorentz5Momentum>::const_iterator p =
	  P.begin() + 2; p != P.end(); ++p )
    m.push_back(p->mass());

  size_t n = P.size() - 2;
  vector<Energy> M(n-1);
  M[0] = Ecm;

  vector<Lorentz5Momentum> Q(n-1);
  Q[0] = Lorentz5Momentum(M[0]);

  for ( size_t i = 2; i <= n-1; ++i ) {
    for ( size_t k = i; k <= n; ++k )
      Q[i-1] += P[k+1];
    M[i-1] = Q[i-1].m();
  }

  double weight = invertIntermediates(M,m,r);

  for ( size_t i = 2; i <= n; ++i ) {
    Lorentz5Momentum p = P[i];
    p.boost(-Q[i-2].boostVector());
    r[n-6+2*i] = (p.cosTheta()+1.)/2.;
    double phi = p.phi();
    if ( phi < 0. )
      phi = 2.*Constants::pi + phi;
    r[n-5+2*i] = phi/(2.*Constants::pi);
  }

  return weight;

}


double FlatInvertiblePhasespace::generateTwoToNKinematics(const double* r,
							  vector<Lorentz5Momentum>& momenta) {

  double weight = generateKinematics(momenta,sqrt(lastXCombPtr()->lastSHat()),r);

  fillDiagramWeights();

  return weight;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FlatInvertiblePhasespace::persistentOutput(PersistentOStream &) const {}

void FlatInvertiblePhasespace::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FlatInvertiblePhasespace,MatchboxPhasespace>
  describeHerwigFlatInvertiblePhasespace("Herwig::FlatInvertiblePhasespace", "Herwig.so");

void FlatInvertiblePhasespace::Init() {

  static ClassDocumentation<FlatInvertiblePhasespace> documentation
    ("FlatInvertiblePhasespace implements flat, invertible phase space generation.");

}

