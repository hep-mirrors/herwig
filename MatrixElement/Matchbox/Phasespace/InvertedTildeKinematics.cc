// -*- C++ -*-
//
// InvertedTildeKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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

#include "Herwig++/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

using namespace Herwig;

InvertedTildeKinematics::InvertedTildeKinematics() 
  : HandlerBase(), theJacobian(0.0), thePtCut(0.0*GeV) {}

InvertedTildeKinematics::~InvertedTildeKinematics() {}

void InvertedTildeKinematics::dumpInfo(const string& prefix) const {
  generator()->log() << prefix << fullName()
		     << " [" << this << "]\n";
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

Lorentz5Momentum InvertedTildeKinematics::getKt (const Lorentz5Momentum& p1,
						 const Lorentz5Momentum& p2,
						 Energy pt,
						 double phi) const {

  Boost beta = (p1+p2).findBoostToCM();
      
  Lorentz5Momentum p1c = p1;

  if (beta.mag2() > Constants::epsilon) {
    p1c.boost(beta);
  }
      
  Lorentz5Momentum k (0.*GeV,0.*GeV,0.*GeV,0.*GeV);
      
  double ct = p1c.vect().unit().z();
  double st = sqrt(1.-ct*ct);
      
  double cphi = cos(phi);
  double sphi = sqrt(1.-cphi*cphi);
  if (phi  > Constants::pi) sphi = -sphi;
      
  if (st > Constants::epsilon) {
    double cchi = p1c.vect().unit().x()/st;
    double schi = p1c.vect().unit().y()/st;
    k.setX((cphi*cchi*ct-sphi*schi)*pt);
    k.setY((cphi*schi*ct+sphi*cchi)*pt);
    k.setZ(-cphi*st*pt);
  } else {
    k.setX(pt*cphi);
    k.setY(pt*sphi);
    k.setZ(0.*GeV);
  }
      
  if (beta.mag2() > Constants::epsilon)
    k.boost(-beta);
      
  return k;

}

Energy InvertedTildeKinematics::lastScale() const {
  if ( ( theDipole->bornEmitter() < 2 && theDipole->bornSpectator() > 1 ) ||
       ( theDipole->bornEmitter() > 1 && theDipole->bornSpectator() < 2 ) ) {
    return -(bornEmitterMomentum()-bornSpectatorMomentum()).m();
  }
  return (bornEmitterMomentum()+bornSpectatorMomentum()).m();
}

pair<Energy,double> InvertedTildeKinematics::generatePtZ(double& jac, const double * r) const {

  double kappaMin = 
    ptCut() != ZERO ?
    sqr(ptCut()/ptMax()) :
    sqr(0.1*GeV/GeV);

  double kappa;

  using namespace RandomHelpers;

  if ( ptCut() > ZERO ) {
    pair<double,double> kw =
      generate(inverse(0.,kappaMin,1.),r[0]);
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

  pair<double,double> zw =
    generate(inverse(0.,zLims.first,zLims.second)+
	     inverse(1.,zLims.first,zLims.second),r[1]);

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
describeInvertedTildeKinematics("Herwig::InvertedTildeKinematics", "HwMatchbox.so");
