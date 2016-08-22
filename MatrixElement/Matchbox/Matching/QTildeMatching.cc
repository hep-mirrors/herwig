// -*- C++ -*-
//
// QTildeMatching.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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

#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/TildeKinematics.h"

using namespace Herwig;

QTildeMatching::QTildeMatching() 
  : theCorrectForXZMismatch(true) {}

QTildeMatching::~QTildeMatching() {}

IBPtr QTildeMatching::clone() const {
  return new_ptr(*this);
}

IBPtr QTildeMatching::fullclone() const {
  return new_ptr(*this);
}

void QTildeMatching::checkCutoff() {
  if ( showerTildeKinematics() ) {
    showerTildeKinematics()->
      prepare(realCXComb(),bornCXComb());
    showerTildeKinematics()->dipole(dipole());
    showerTildeKinematics()->getShowerVariables();
  }
}

void QTildeMatching::getShowerVariables() {

  // already filled from checkCutoff in this case
  if ( showerTildeKinematics() )
    return;

  // get the shower variables
  calculateShowerVariables();

  // check for the cutoff
  dipole()->isAboveCutoff(isAboveCutoff());

  // get the hard scale
  dipole()->showerHardScale(hardScale());

  // check for phase space
  dipole()->isInShowerPhasespace(isInShowerPhasespace());

}

bool QTildeMatching::isInShowerPhasespace() const {

  assert((theQTildeSudakov->cutOffOption() == 0 || theQTildeSudakov->cutOffOption() == 2) && 
	 "implementation only provided for default and pt cutoff");

  Energy qtildeHard = ZERO;

  Energy qtilde = dipole()->showerScale();
  assert(!dipole()->showerParameters().empty());
  double z = dipole()->showerParameters()[0];

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

  if ( pt2 < max(theQTildeSudakov->pT2min(),sqr(safeCut()) ))
    return false;

  bool hardVeto = restrictPhasespace() && sqrt(pt2) >= dipole()->showerHardScale();
  return qtilde <= qtildeHard && !hardVeto;

}

bool QTildeMatching::isAboveCutoff() const {

  assert((theQTildeSudakov->cutOffOption() == 0 || theQTildeSudakov->cutOffOption() == 2) && 
	 "implementation only provided for default and pt cutoff");
  Energy qtilde = dipole()->showerScale();
  assert(!dipole()->showerParameters().empty());
  double z = dipole()->showerParameters()[0];
  Energy Qg = theQTildeSudakov->kinScale();
  if ( dipole()->bornEmitter() > 1 ) {
    Energy mu = max(Qg,realCXComb()->meMomenta()[dipole()->realEmitter()].mass());
    if ( bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id() == ParticleID::g )
      return sqr(z*(1.-z)*qtilde) - sqr(mu) >= 
             max(theQTildeSudakov->pT2min(),sqr(safeCut()));
    else
      return sqr(z*(1.-z)*qtilde) - sqr((1.-z)*mu) - z*sqr(Qg) >= 
             max(theQTildeSudakov->pT2min(),sqr(safeCut()));
  }
  if ( dipole()->bornEmitter() < 2 ) {
    return
      sqr((1.-z)*qtilde) - z*sqr(Qg) >= 
      max(theQTildeSudakov->pT2min(),sqr(safeCut()));
  }
  return false;
}

CrossSection QTildeMatching::dSigHatDR() const {

  assert(!dipole()->showerParameters().empty());
  pair<Energy2,double> vars = 
    make_pair(sqr(dipole()->showerScale()),
	      dipole()->showerParameters()[0]);

  pair<int,int> ij(dipole()->bornEmitter(),
		   dipole()->bornSpectator());
  double ccme2 =
     dipole()->underlyingBornME()->largeNColourCorrelatedME2(ij,theLargeNBasis);
   
   if(ccme2==0.)return 0.*nanobarn;
      
   double lnme2=dipole()->underlyingBornME()->largeNME2(theLargeNBasis);
   if(lnme2==0){
     generator()->log() <<"\nQTildeMatching: ";
     generator()->log() <<"\n  largeNME2 is ZERO, while largeNColourCorrelatedME2 is not ZERO." ;
     generator()->log() <<"\n  This is too seriuos.\n" ;
     generator()->log() << Exception::runerror;
   }
  
   ccme2 *=
     dipole()->underlyingBornME()->me2() /lnme2;
     

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
  
  xme2 *= dipole()->realEmissionME()->finalStateSymmetry() /
    dipole()->underlyingBornME()->finalStateSymmetry(); 

  // take care of mismatch between z and x as we are approaching the
  // hard phase space boundary
  // TODO get rid of this useless scale option business and simplify PDF handling in here
  if ( dipole()->bornEmitter() < 2 && theCorrectForXZMismatch ) {
    Energy2 emissionScale = ZERO;
    if ( emissionScaleInSubtraction() == showerScale ) {
      emissionScale = showerFactorizationScale();
    } else if ( emissionScaleInSubtraction() == realScale ) {
      emissionScale = dipole()->realEmissionME()->lastScale();
    } else if ( emissionScaleInSubtraction() == bornScale ) {
      emissionScale = dipole()->underlyingBornME()->lastScale();
    }
    double xzMismatch = 
      dipole()->subtractionParameters()[0] / dipole()->showerParameters()[0];
    double realCorrectedPDF = 
      dipole()->bornEmitter() == 0 ?
      dipole()->realEmissionME()->pdf1(emissionScale,theExtrapolationX,
				       xzMismatch) :
      dipole()->realEmissionME()->pdf2(emissionScale,theExtrapolationX,
				       xzMismatch);
    double realPDF = 
      dipole()->bornEmitter() == 0 ?
      dipole()->realEmissionME()->pdf1(emissionScale,theExtrapolationX,1.0) :
      dipole()->realEmissionME()->pdf2(emissionScale,theExtrapolationX,1.0);
    if ( realPDF == 0.0 || realCorrectedPDF == 0.0 )
      return ZERO;
    xme2 *= realCorrectedPDF / realPDF;
  }

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

  if ( profileScales() )
    xme2 *= profileScales()->hardScaleProfile(dipole()->showerHardScale(),sqrt(pt2));

  CrossSection res = 
    sqr(hbarc) * 
    realXComb()->jacobian() * 
    subtractionScaleWeight() *
    xme2 /
    (2. * realXComb()->lastSHat());

  return res;

}

double QTildeMatching::me2() const {
  throw Exception() << "QTildeMatching::me2(): Not intented to use. Disable the ShowerApproximationGenerator."
		    << Exception::runerror;
  return 0.;
}

void QTildeMatching::calculateShowerVariables() const {

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

  // the light-cone condition is numerically not very stable, so we
  // explicitly push it on the light-cone here
  n.setMass(ZERO);
  n.rescaleEnergy();

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
  // allow small violations (numerical inaccuracies)
  if ( z <= 0 && z >= -1e-6 ) {
    z = std::numeric_limits<double>::epsilon();
  } else if ( z >= 1 && z <= 1+1e-6 ) {
    z = 1-std::numeric_limits<double>::epsilon();
  }

  Energy2 qtilde2 = ZERO;
  Energy2 q2 = ZERO;

  if ( dipole()->bornEmitter() > 1 ) {
    q2 = 
      (realCXComb()->meMomenta()[dipole()->realEmitter()] + realCXComb()->meMomenta()[dipole()->realEmission()]).m2();
    qtilde2 = (q2 - bornCXComb()->meMomenta()[dipole()->bornEmitter()].m2())/(z*(1.-z));
  } else {
    q2 = 
      -(realCXComb()->meMomenta()[dipole()->realEmitter()] - realCXComb()->meMomenta()[dipole()->realEmission()]).m2();
    qtilde2 = (q2 + bornCXComb()->meMomenta()[dipole()->bornEmitter()].m2())/(1.-z);
  }
  if ( qtilde2 < ZERO ) {
    qtilde2 = ZERO;
  }

  assert(qtilde2 >= ZERO && z > 0.0 && z < 1.0);

  dipole()->showerScale(sqrt(qtilde2));
  dipole()->showerParameters().resize(1);
  dipole()->showerParameters()[0] = z;

}

double QTildeMatching::splitFn(const pair<Energy2,double>& vars) const {
  const Energy2& qtilde2 = vars.first;
  const double& z = vars.second;
  double Nc = SM().Nc();

  // final state branching
  if ( dipole()->bornEmitter() > 1 ) {
    // final state quark quark branching
    if ( abs(bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id()) < 7 ) {
      Energy m = bornCXComb()->mePartonData()[dipole()->bornEmitter()]->hardProcessMass();
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
	Energy m = realCXComb()->mePartonData()[dipole()->realEmission()]->hardProcessMass();
	return (1./2.)*(1.-2.*z*(1.-z)+2.*sqr(m)/(z*(1.-z)*qtilde2));
      }
    }
    // final state squark branching
    if ((abs(bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id()) > 1000000 && 
	 abs(bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id()) < 1000007) ||
	(abs(bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id()) > 2000000 && 
	 abs(bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id()) < 2000007)){
      Energy m = bornCXComb()->mePartonData()[dipole()->bornEmitter()]->hardProcessMass();
      return ((sqr(Nc)-1.)/Nc)*(z-sqr(m)/(z*qtilde2))/(1.-z);
    }
    // final state gluino branching
    if (bornCXComb()->mePartonData()[dipole()->bornEmitter()]->id() == 1000021){
      Energy m = bornCXComb()->mePartonData()[dipole()->bornEmitter()]->hardProcessMass();
      return Nc*(1.+sqr(z)-2.*sqr(m)/(z*qtilde2))/(1.-z);
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

void QTildeMatching::doinit() {
  assert(theShowerHandler && theQTildeFinder && theQTildeSudakov);
  theShowerHandler->init();
  theQTildeFinder->init();
  theQTildeSudakov->init();
  hardScaleFactor(theShowerHandler->hardScaleFactor());
  factorizationScaleFactor(theShowerHandler->factorizationScaleFactor());
  renormalizationScaleFactor(theShowerHandler->renormalizationScaleFactor());
  profileScales(theShowerHandler->profileScales());
  restrictPhasespace(theShowerHandler->restrictPhasespace());
  hardScaleIsMuF(theShowerHandler->hardScaleIsMuF());
  ShowerApproximation::doinit();
}

void QTildeMatching::doinitrun() {
  assert(theShowerHandler && theQTildeFinder && theQTildeSudakov);
  theShowerHandler->initrun();
  theQTildeFinder->initrun();
  theQTildeSudakov->initrun();
  ShowerApproximation::doinitrun();
}


void QTildeMatching::persistentOutput(PersistentOStream & os) const {
  os << theQTildeFinder << theQTildeSudakov
     << theShowerHandler << theCorrectForXZMismatch;
}

void QTildeMatching::persistentInput(PersistentIStream & is, int) {
  is >> theQTildeFinder >> theQTildeSudakov
     >> theShowerHandler >> theCorrectForXZMismatch;
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

  static Reference<QTildeMatching,ShowerHandler> interfaceShowerHandler
    ("ShowerHandler",
     "",
     &QTildeMatching::theShowerHandler, false, false, true, true, false);

  static Switch<QTildeMatching,bool> interfaceCorrectForXZMismatch
    ("CorrectForXZMismatch",
     "Correct for x/z mismatch near hard phase space boundary.",
     &QTildeMatching::theCorrectForXZMismatch, true, false, false);
  static SwitchOption interfaceCorrectForXZMismatchYes
    (interfaceCorrectForXZMismatch,
     "Yes",
     "Include the correction factor.",
     true);
  static SwitchOption interfaceCorrectForXZMismatchNo
    (interfaceCorrectForXZMismatch,
     "No",
     "Do not include the correction factor.",
     false);

}

