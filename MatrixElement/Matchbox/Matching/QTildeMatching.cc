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
  if ( !bornCXComb()->mePartonData()[0]->coloured() &&
       !bornCXComb()->mePartonData()[1]->coloured() ) {
    Energy maxPt = (bornCXComb()->meMomenta()[0] + bornCXComb()->meMomenta()[1]).m();
    maxPt *= hardScaleFactor();
    return maxPt;
  }
  Energy maxPt = -1.0*GeV;
  vector<Lorentz5Momentum>::const_iterator p = 
    bornCXComb()->meMomenta().begin() + 2;
  cPDVector::const_iterator pp = 
    bornCXComb()->mePartonData().begin() + 2;
  for ( ; p != bornCXComb()->meMomenta().end(); ++p, ++pp )
    if ( (**pp).coloured() )
      maxPt = max(maxPt,p->mt());
  if ( maxPt < ZERO )
    maxPt = (bornCXComb()->meMomenta()[0] + bornCXComb()->meMomenta()[1]).m();
  maxPt *= hardScaleFactor();
  return maxPt;
}

double QTildeMatching::hardScaleProfile(Energy Q, Energy q) const {
  if ( q <= Q )
    return 1.;
  return 0.;
}

bool QTildeMatching::isInShowerPhasespace() const {

  assert(theQTildeSudakov->cutOffOption() == 0 && "implementation only provided for default cutoff");

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


  Energy Qg = theQTildeSudakov->kinScale();
  Energy2 pt2 = ZERO;
  if ( dipole()->bornEmitter() > 1 ) {
    Energy mu = max(Qg,realCXComb()->meMomenta()[dipole()->realEmitter()].mass());
    if ( bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id() == ParticleID::g )
      pt2 = sqr(z*(1.-z)*qtilde) - sqr(mu);
    else
      pt2 = sqr(z*(1.-z)*qtilde) - sqr((1.-z)*mu) - z*sqr(Qg);
  }
  if ( dipole()->bornEmitter() < 2 ) {
    pt2 = sqr((1.-z)*qtilde) - z*sqr(Qg);
  }

  if ( pt2 < ZERO )
    return false;

  return qtilde <= qtildeHard && sqrt(pt2) < hardScale();

}

bool QTildeMatching::isAboveCutoff() const {
  assert(theQTildeSudakov->cutOffOption() == 0 && "implementation only provided for default cutoff");
  pair<Energy2,double> vars = getShowerVariables();
  Energy qtilde = sqrt(vars.first);
  double z = vars.second;
  Energy Qg = theQTildeSudakov->kinScale();
  if ( dipole()->bornEmitter() > 1 ) {
    Energy mu = max(Qg,realCXComb()->meMomenta()[dipole()->realEmitter()].mass());
    if ( bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id() == ParticleID::g )
      return sqr(z*(1.-z)*qtilde) - sqr(mu) >= ZERO;
    else
      return sqr(z*(1.-z)*qtilde) - sqr((1.-z)*mu) - z*sqr(Qg) >= ZERO;
  }
  if ( dipole()->bornEmitter() < 2 ) {
    return
      sqr((1.-z)*qtilde) - z*sqr(Qg) >= ZERO;
  }
  return false;
}

CrossSection QTildeMatching::dSigHatDR() const {

  pair<Energy2,double> vars = getShowerVariables();

  pair<int,int> ij(dipole()->bornEmitter(),
		   dipole()->bornSpectator());
  double ccme2 = 
    theShowerKernels ?
    dipole()->underlyingBornME()->largeNColourCorrelatedME2(ij,theLargeNBasis) :
    dipole()->underlyingBornME()->colourCorrelatedME2(ij);

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

  Energy qtilde = sqrt(vars.first);
  double z = vars.second;
  Energy2 pt2 = ZERO;
  Energy Qg = theQTildeSudakov->kinScale();
  if ( dipole()->bornEmitter() > 1 ) {
    Energy mu = max(Qg,realCXComb()->meMomenta()[dipole()->realEmitter()].mass());
    if ( bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id() == ParticleID::g )
      pt2 = sqr(z*(1.-z)*qtilde) - sqr(mu);
    else
      pt2 = sqr(z*(1.-z)*qtilde) - sqr((1.-z)*mu) - z*sqr(Qg);
  }
  if ( dipole()->bornEmitter() < 2 ) {
    pt2 = sqr((1.-z)*qtilde) - z*sqr(Qg);
  }
  assert(pt2 >= ZERO);
  xme2 *= hardScaleProfile(hardScale(),sqrt(pt2));

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

  Lorentz5Momentum n;
  Energy2 Q2 = ZERO;

  const Lorentz5Momentum& pb = bornCXComb()->meMomenta()[dipole()->bornEmitter()];
  const Lorentz5Momentum& pc = bornCXComb()->meMomenta()[dipole()->bornSpectator()];

  if ( dipole()->bornEmitter() > 1 ) {
    Q2 = (pb+pc).m2();
  } else {
    Q2 = -(pb-pc).m2();
  }

  if ( dipole()->bornEmitter() > 1 && dipole()->bornSpectator() > 1 ) {
    double b = sqr(bornCXComb()->meMomenta()[dipole()->bornEmitter()].m())/Q2;
    double c = sqr(bornCXComb()->meMomenta()[dipole()->bornSpectator()].m())/Q2;
    double lambda = sqrt(1.+sqr(b)+sqr(c)-2.*b-2.*c-2.*b*c);
    n = (1.-0.5*(1.-b+c-lambda))*pc - 0.5*(1.-b+c-lambda)*pb;
  }

  if ( dipole()->bornEmitter() > 1 && dipole()->bornSpectator() < 2 ) {
    n = bornCXComb()->meMomenta()[dipole()->bornSpectator()];
  }

  if ( dipole()->bornEmitter() < 2 && dipole()->bornSpectator() > 1 ) {
    double c = sqr(bornCXComb()->meMomenta()[dipole()->bornSpectator()].m())/Q2;
    n = (1.+c)*pc - c*pb;
  }

  if ( dipole()->bornEmitter() < 2 && dipole()->bornSpectator() < 2 ) {
    n = bornCXComb()->meMomenta()[dipole()->bornSpectator()];
  }

  assert(abs(n.m()/GeV) < 1e-9);

  double z = 0.0;

  if ( dipole()->bornEmitter() > 1 ) {
    z = 1. - 
    (n*realCXComb()->meMomenta()[dipole()->realEmission()])/
    (n*bornCXComb()->meMomenta()[dipole()->bornEmitter()]);
  } else {
    z = 1. - 
    (n*realCXComb()->meMomenta()[dipole()->realEmission()])/
    (n*realCXComb()->meMomenta()[dipole()->realEmitter()]);
  }

  Energy2 qtilde2 = ZERO;

  if ( dipole()->bornEmitter() > 1 ) {
    qtilde2 = (Q2 - bornCXComb()->meMomenta()[dipole()->bornEmitter()].m2())/(z*(1.-z));
  } else {
    qtilde2 = (Q2 + bornCXComb()->meMomenta()[dipole()->bornEmitter()].m2())/(1.-z);
  }

  assert(qtilde2 >= ZERO && z >= 0.0 && z <= 1.0);

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
  os << theQTildeFinder << theQTildeSudakov;
}

void QTildeMatching::persistentInput(PersistentIStream & is, int) {
  is >> theQTildeFinder >> theQTildeSudakov;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<QTildeMatching,Herwig::ShowerApproximation>
  describeHerwigQTildeMatching("Herwig::QTildeMatching", "HwShower.so HwQTildeMatching.so");

void QTildeMatching::Init() {

  static ClassDocumentation<QTildeMatching> documentation
    ("QTildeMatching implements NLO matching with the default shower.");

  static Reference<QTildeMatching,QTildeFinder> interfaceQTildeFinder
    ("QTildeFinder",
     "Set the partner finder to calculate hard scales.",
     &QTildeMatching::theQTildeFinder, false, false, true, false, false);

  static Reference<QTildeMatching,QTildeSudakov> interfaceQTildeSudakov
    ("QTildeSudakov",
     "Set the partner finder to calculate hard scales.",
     &QTildeMatching::theQTildeSudakov, false, false, true, false, false);

}

