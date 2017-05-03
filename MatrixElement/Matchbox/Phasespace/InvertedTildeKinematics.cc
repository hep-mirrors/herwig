// -*- C++ -*-
//
// InvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the InvertedTildeKinematics class.
//

#include <limits>

#include "InvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Rebinder.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

using namespace Herwig;

InvertedTildeKinematics::InvertedTildeKinematics() 
  : HandlerBase(), theJacobian(0.0), thePtCut(0.0*GeV) {}

InvertedTildeKinematics::~InvertedTildeKinematics() {}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

Lorentz5Momentum InvertedTildeKinematics::getKt(const Lorentz5Momentum& p1,
						const Lorentz5Momentum& p2,
						Energy pt,
						double phi,
						bool spacelike) const {

  Lorentz5Momentum P;
  if ( !spacelike )
    P = p1 + p2;
  else
    P = p1 - p2;

  Energy2 Q2 = abs(P.m2());

  Lorentz5Momentum Q = 
    !spacelike ? 
    Lorentz5Momentum(ZERO,ZERO,ZERO,sqrt(Q2),sqrt(Q2)) :
    Lorentz5Momentum(ZERO,ZERO,sqrt(Q2),ZERO,-sqrt(Q2));

  if ( spacelike && Q.z() < P.z() )
    Q.setZ(-Q.z());

  bool boost =
    abs((P-Q).vect().mag2()/GeV2) > 1e-10 ||
    abs((P-Q).t()/GeV) > 1e-5;
  boost &= (P*Q-Q.mass2())/GeV2 > 1e-8;

  Lorentz5Momentum inFrame1;
  if ( boost )
    inFrame1 = p1 + ((P*p1-Q*p1)/(P*Q-Q.mass2()))*(P-Q);
  else
    inFrame1 = p1;

  Energy ptx = inFrame1.x();
  Energy pty = inFrame1.y();
  Energy q = 2.*inFrame1.z();

  Energy Qp = sqrt(4.*(sqr(ptx)+sqr(pty))+sqr(q));
  Energy Qy = sqrt(4.*sqr(pty)+sqr(q));

  double cPhi = cos(phi);
  double sPhi = sqrt(1.-sqr(cPhi));
  if ( phi > Constants::pi )
    sPhi = -sPhi;

  Lorentz5Momentum kt;

  if ( !spacelike ) {
    kt.setT(ZERO);
    kt.setX(pt*Qy*cPhi/Qp);
    kt.setY(-pt*(4*ptx*pty*cPhi/Qp+q*sPhi)/Qy);
    kt.setZ(2.*pt*(-ptx*q*cPhi/Qp + pty*sPhi)/Qy);
  } else {
    kt.setT(2.*pt*(ptx*q*cPhi+pty*Qp*sPhi)/(q*Qy));
    kt.setX(pt*(Qp*q*cPhi+4.*ptx*pty*sPhi)/(q*Qy));
    kt.setY(pt*Qy*sPhi/q);
    kt.setZ(ZERO);
  }

  if ( boost )
    kt = kt + ((P*kt-Q*kt)/(P*Q-Q.mass2()))*(P-Q);
  kt.setMass(-pt);
  kt.rescaleRho();

  return kt;

}

Energy InvertedTildeKinematics::lastScale() const {
  if ( ( theDipole->bornEmitter() < 2 && theDipole->bornSpectator() > 1 ) ||
       ( theDipole->bornEmitter() > 1 && theDipole->bornSpectator() < 2 ) ) {
    return -(bornEmitterMomentum()-bornSpectatorMomentum()).m();
  }
  return (bornEmitterMomentum()+bornSpectatorMomentum()).m();
}

pair<Energy,double> InvertedTildeKinematics::generatePtZ(double& jac, const double * r,
							 double pow, vector<double>* ) const {

  double kappaMin = 
    ptCut() != ZERO ?
    sqr(ptCut()/ptMax()) :
    sqr(0.1*GeV/GeV);

  double kappa;

  using namespace RandomHelpers;

  if ( ptCut() > ZERO ) {
    pair<double,double> kw = pow==1. ?
      generate(inverse(0.,kappaMin,1.),r[0]) :
      generate(power(0.,-pow,kappaMin,1.),r[0]);
    kappa = kw.first;
    jac *= kw.second;
  } else {
    pair<double,double> kw =
      generate((piecewise(),
		flat(1e-4,kappaMin),
		match(inverse(0.,kappaMin,1.))),r[0]);
    kappa = kw.first;
    jac *= kw.second;
  }

  Energy pt = sqrt(kappa)*ptMax();

  pair<double,double> zLims = zBounds(pt);

  pair<double,double> zw(0,0);// =
  //  generate(inverse(0.,zLims.first,zLims.second)+
	//     inverse(1.,zLims.first,zLims.second),r[1]);

  // FlatZ = 1
  if ( theDipole->samplingZ() == 1 ) {
    zw = generate(flat(zLims.first,zLims.second),r[1]);
  }
  // OneOverZ = 2
  if ( theDipole->samplingZ() == 2 ) {
    zw = generate(inverse(0.0,zLims.first,zLims.second),r[1]);
  }
  // OneOverOneMinusZ = 3
  if ( theDipole->samplingZ() == 3 ) {
    zw = generate(inverse(1.0,zLims.first,zLims.second),r[1]);
  }
  // OneOverZOneMinusZ = 4
  if ( theDipole->samplingZ() == 4 ) {
    zw = generate(inverse(0.0,zLims.first,zLims.second) +
                                inverse(1.0,zLims.first,zLims.second),r[1]);
  }

  double z = zw.first;
  jac *= zw.second;

  jac *= sqr(ptMax()/lastScale());

  return make_pair(pt,z);

}

void InvertedTildeKinematics::rebind(const TranslationMap & trans) {
  theDipole = trans.translate(theDipole);
  HandlerBase::rebind(trans);
}

IVector InvertedTildeKinematics::getReferences() {
  IVector ret = HandlerBase::getReferences();
  ret.push_back(theDipole);
  return ret;
}

void InvertedTildeKinematics::persistentOutput(PersistentOStream & os) const {
  os << theDipole << theRealXComb << theBornXComb
     << ounit(theRealEmitterMomentum,GeV) << ounit(theRealEmissionMomentum,GeV)
     << ounit(theRealSpectatorMomentum,GeV) << theJacobian
     << ounit(thePtCut,GeV);
}

void InvertedTildeKinematics::persistentInput(PersistentIStream & is, int) {
  is >> theDipole >> theRealXComb >> theBornXComb
     >> iunit(theRealEmitterMomentum,GeV) >> iunit(theRealEmissionMomentum,GeV)
     >> iunit(theRealSpectatorMomentum,GeV) >> theJacobian
     >> iunit(thePtCut,GeV);
}

void InvertedTildeKinematics::Init() {

  static ClassDocumentation<InvertedTildeKinematics> documentation
    ("InvertedTildeKinematics is the base class for the inverted 'tilde' "
     "kinematics being used for subtraction terms in the "
     "formalism of Catani and Seymour.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<InvertedTildeKinematics,HandlerBase>
describeInvertedTildeKinematics("Herwig::InvertedTildeKinematics", "Herwig.so");
