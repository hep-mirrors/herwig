// -*- C++ -*-
//
// QTildeMatching.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeMatching class.
//

#include "QTildeMatching.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

using namespace Herwig;

QTildeMatching::QTildeMatching() {}

QTildeMatching::~QTildeMatching() {}

IBPtr QTildeMatching::clone() const {
  return new_ptr(*this);
}

IBPtr QTildeMatching::fullclone() const {
  return new_ptr(*this);
}

Energy QTildeMatching::hardScale() const {
  // should be the same as for the dipole shower but needs checking;
  // NB this is the pt veto scale as the `hard' scale is anyway fixed
  return ShowerApproximation::hardScale();
}

double QTildeMatching::hardScaleProfile(Energy, Energy) const {
  // this needs support in the shower - leave as is for the time being
  return 1.;
}

bool QTildeMatching::isInShowerPhasespace() const {

  if ( !isAboveCutoff() )
    return false;

  if ( dipole()->lastPt() > hardScale() )
    return false;

  Energy qtildeHard = ZERO;

  pair<Energy2,double> vars = getShowerVariables();
  Energy qtilde = sqrt(vars.first);
  double z = vars.second;

  // FF
  if ( dipole()->bornEmitter() > 1 && dipole()->bornSpectator() > 1 ) {
    qtildeHard = 
      theQTildeFinder->
      calculateFinalFinalScales(bornCXComb()->meMomenta()[dipole()->bornEmitter()],
				bornCXComb()->meMomenta()[dipole()->bornSpectator()],
				bornCXComb()->mePartonData()[dipole()->bornEmitter()]->iColour() == PDT::Colour3).first;
  }

  // FI
  if ( dipole()->bornEmitter() > 1 && dipole()->bornSpectator() < 2 ) {
    qtildeHard = 
      theQTildeFinder->
      calculateInitialFinalScales(bornCXComb()->meMomenta()[dipole()->bornSpectator()],
				  bornCXComb()->meMomenta()[dipole()->bornEmitter()],false).second;
  }

  // IF
  if ( dipole()->bornEmitter() < 2 && dipole()->bornSpectator() > 1 ) {
    qtildeHard = 
      theQTildeFinder->
      calculateInitialFinalScales(bornCXComb()->meMomenta()[dipole()->bornEmitter()],
				  bornCXComb()->meMomenta()[dipole()->bornSpectator()],false).first;
    if ( z < (dipole()->bornEmitter() == 0 ? bornCXComb()->lastX1() : bornCXComb()->lastX2()) )
      return false;
  }

  // II
  if ( dipole()->bornEmitter() < 2 && dipole()->bornSpectator() < 2 ) {
    qtildeHard = 
      theQTildeFinder->
      calculateInitialInitialScales(bornCXComb()->meMomenta()[dipole()->bornEmitter()],
				    bornCXComb()->meMomenta()[dipole()->bornSpectator()]).first;
    if ( z < (dipole()->bornEmitter() == 0 ? bornCXComb()->lastX1() : bornCXComb()->lastX2()) )
      return false;
  }

  return qtilde <= qtildeHard;

}

bool QTildeMatching::isAboveCutoff() const {
  assert(theQTildeSudakov->cutOffOption() == 0 && "implementation only provided for default cutoff");
  pair<Energy2,double> vars = getShowerVariables();
  Energy qtilde = sqrt(vars.first);
  double z = vars.second;
  Energy Qg = theQTildeSudakov->kinScale();
  if ( dipole()->bornEmitter() > 1 ) {
    Energy mu = max(Qg,realCXComb()->meMomenta()[dipole()->realEmitter()].mass());
    if ( abs(realCXComb()->mePartonData()[dipole()->realEmission()]->id()) < 7 &&
	 bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id() == ParticleID::g )
      mu = realCXComb()->meMomenta()[dipole()->realEmitter()].mass();
    return sqr(z*(1.-z)*qtilde) >= sqr((1.-z)*mu)+z*sqr(Qg);
  }
  if ( dipole()->bornEmitter() < 2 ) {
    return
      z <= 1.+Qg/(2.*qtilde) - sqrt(sqr(1.+Qg/(2.*qtilde))-1.);
  }
  return false;
}

CrossSection QTildeMatching::dSigHatDR() const {

  pair<Energy2,double> vars = getShowerVariables();

  pair<int,int> ij(dipole()->bornEmitter(),
		   dipole()->bornSpectator());
  double ccme2 = 
    dipole()->underlyingBornME()->largeNColourCorrelatedME2(ij,theLargeNBasis);

  Energy2 prop = ZERO;
  if ( dipole()->bornEmitter() > 1 ) {
    prop =
      (realCXComb()->meMomenta()[dipole()->realEmitter()] +
       realCXComb()->meMomenta()[dipole()->realEmission()]).m2()
      - bornCXComb()->meMomenta()[dipole()->bornEmitter()].m2();
  } else {
    prop = 
      2.*vars.second*(realCXComb()->meMomenta()[dipole()->realEmitter()]*
		      realCXComb()->meMomenta()[dipole()->realEmission()]);
  }

  // note alphas included downstream from subtractionScaleWeight()
  double xme2 = -8.*Constants::pi*ccme2*splitFn(vars)*realXComb()->lastSHat()/prop;
  xme2 *= 
    pow(realCXComb()->lastSHat() / bornCXComb()->lastSHat(),
	bornCXComb()->mePartonData().size()-4.);

  double bornPDF = bornPDFWeight(dipole()->underlyingBornME()->lastScale());
  if ( bornPDF == 0.0 )
    return ZERO;

  xme2 *= bornPDF;

  xme2 *= hardScaleProfile(hardScale(),dipole()->lastPt());

  CrossSection res = 
    sqr(hbarc) * 
    realXComb()->jacobian() * 
    subtractionScaleWeight() *
    xme2 /
    (2. * realXComb()->lastSHat());

  return res;

}

double QTildeMatching::me2() const {
  throw Exception() << "Not intented to use. Disable the ShowerApproximationGenerator."
		    << Exception::abortnow;
  return 0.;
}

pair<Energy2,double> QTildeMatching::getShowerVariables() const {

  Energy2 qtilde2 = ZERO;
  double z = 1. -
    (bornCXComb()->meMomenta()[dipole()->bornSpectator()]*
     realCXComb()->meMomenta()[dipole()->realEmission()]) / 
    (bornCXComb()->meMomenta()[dipole()->bornSpectator()]*
     bornCXComb()->meMomenta()[dipole()->bornEmitter()]);

  // final state branching
  if ( dipole()->bornEmitter() > 1 ) {
    // final state quark quark branching
    if ( abs(bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id()) < 7 ) {
      qtilde2 = 
	sqr(dipole()->lastPt())/sqr(z*(1.-z)) +
	sqr(bornCXComb()->mePartonData()[dipole()->bornEmitter()]->mass())/sqr(z);
    }
    // final state gluon branching
    if ( bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id() == ParticleID::g ) {
      qtilde2 = 
	(sqr(dipole()->lastPt()) +
	 sqr(realCXComb()->mePartonData()[dipole()->realEmitter()]->mass()))/sqr(z*(1.-z));
    }
  }

  // initial state branching
  if ( dipole()->bornEmitter() < 2 ) {
    qtilde2 =
      sqr(dipole()->lastPt())/sqr(z*(1.-z));
  }

  return make_pair(qtilde2,z);

}

double QTildeMatching::splitFn(const pair<Energy2,double>& vars) const {
  const Energy2& qtilde2 = vars.first;
  const double& z = vars.second;
  double Nc = SM().Nc();

  // final state branching
  if ( dipole()->bornEmitter() > 1 ) {
    // final state quark quark branching
    if ( abs(bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id()) < 7 ) {
      Energy m = bornCXComb()->mePartonData()[dipole()->bornEmitter()]->mass();
      return
	((sqr(Nc)-1.)/(2.*Nc))*(1+sqr(z)-2.*sqr(m)/(z*qtilde2))/(1.-z);
    }
    // final state gluon branching
    if ( bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id() == ParticleID::g ) {
      if ( realCXComb()->mePartonData()[dipole()->realEmission()]->id() == ParticleID::g ) {
	// ATTENTION the factor 2 here is intentional as it cancels to the 1/2
	// stemming from the large-N colour correlator
	return 2.*Nc*(z/(1.-z)+(1.-z)/z+z*(1.-z));
      }
      if ( abs(realCXComb()->mePartonData()[dipole()->realEmission()]->id()) < 7 ) {
	Energy m = realCXComb()->mePartonData()[dipole()->realEmission()]->mass();
	return (1./2.)*(1.-2.*z*(1.-z)+2.*sqr(m)/(z*(1.-z)*qtilde2));
      }
    }
  }

  // initial state branching
  if ( dipole()->bornEmitter() < 2 ) {
    // g/g
    if ( realCXComb()->mePartonData()[dipole()->realEmitter()]->id() == ParticleID::g &&
	 realCXComb()->mePartonData()[dipole()->realEmission()]->id() == ParticleID::g ) {
      // see above for factor of 2
      return 2.*Nc*(z/(1.-z)+(1.-z)/z+z*(1.-z));
    }
    // q/q
    if ( abs(realCXComb()->mePartonData()[dipole()->realEmitter()]->id()) < 7 &&
	 realCXComb()->mePartonData()[dipole()->realEmission()]->id() == ParticleID::g ) {
      return
	((sqr(Nc)-1.)/(2.*Nc))*(1+sqr(z))/(1.-z);
    }
    // g/q
    if ( realCXComb()->mePartonData()[dipole()->realEmitter()]->id() == ParticleID::g &&
	 abs(realCXComb()->mePartonData()[dipole()->realEmission()]->id()) < 7 ) {
      return (1./2.)*(1.-2.*z*(1.-z));
    }
    // q/g
    if ( abs(realCXComb()->mePartonData()[dipole()->realEmitter()]->id()) < 7 &&
	 abs(realCXComb()->mePartonData()[dipole()->realEmission()]->id()) < 7 ) {
      return
	((sqr(Nc)-1.)/(2.*Nc))*(1+sqr(1.-z))/z;
    }
  }

  return 0.0;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void QTildeMatching::persistentOutput(PersistentOStream & os) const {
  os << theLargeNBasis << theQTildeFinder << theQTildeSudakov;
}

void QTildeMatching::persistentInput(PersistentIStream & is, int) {
  is >> theLargeNBasis >> theQTildeFinder >> theQTildeSudakov;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<QTildeMatching,Herwig::ShowerApproximation>
  describeHerwigQTildeMatching("Herwig::QTildeMatching", "HwQTildeMatching.so HwShower.so");

void QTildeMatching::Init() {

  static ClassDocumentation<QTildeMatching> documentation
    ("QTildeMatching implements NLO matching with the default shower.");

  static Reference<QTildeMatching,ColourBasis> interfaceLargeNBasis
    ("LargeNBasis",
     "Set the large-N colour basis implementation.",
     &QTildeMatching::theLargeNBasis, false, false, true, true, false);

  static Reference<QTildeMatching,QTildeFinder> interfaceQTildeFinder
    ("QTildeFinder",
     "Set the partner finder to calculate hard scales.",
     &QTildeMatching::theQTildeFinder, false, false, true, false, false);

  static Reference<QTildeMatching,QTildeSudakov> interfaceQTildeSudakov
    ("QTildeSudakov",
     "Set the partner finder to calculate hard scales.",
     &QTildeMatching::theQTildeSudakov, false, false, true, false, false);

}

