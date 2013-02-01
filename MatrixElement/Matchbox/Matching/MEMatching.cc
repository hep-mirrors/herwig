// -*- C++ -*-
//
// MEMatching.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEMatching class.
//

#include "MEMatching.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/PDT/EnumParticles.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

using namespace Herwig;

MEMatching::MEMatching() {}

MEMatching::~MEMatching() {}

IBPtr MEMatching::clone() const {
  return new_ptr(*this);
}

IBPtr MEMatching::fullclone() const {
  return new_ptr(*this);
}

double MEMatching::channelWeight(int emitter, int emission, int spectator) const {
  // do the most simple thing for the time being; needs fixing later
  if ( realCXComb()->mePartonData()[emission]->id() == ParticleID::g ) {
    Energy2 pipk = 
      realCXComb()->meMomenta()[emitter] * realCXComb()->meMomenta()[spectator];
    Energy2 pipj = 
      realCXComb()->meMomenta()[emitter] * realCXComb()->meMomenta()[emission];
    Energy2 pjpk = 
      realCXComb()->meMomenta()[emission] * realCXComb()->meMomenta()[spectator];
    return GeV2 * pipk / ( pipj * ( pipj + pjpk ) );
  }
  return
    GeV2 / (realCXComb()->meMomenta()[emitter] * realCXComb()->meMomenta()[emission]);
}

double MEMatching::channelWeight() const {
  double currentChannel = channelWeight(dipole()->realEmitter(),
					dipole()->realEmission(),
					dipole()->realSpectator());
  if ( currentChannel == 0. )
    return 0.;
  double sum = 0.;
  for ( vector<Ptr<SubtractionDipole>::ptr>::const_iterator dip =
	  dipole()->partnerDipoles().begin();
	dip != dipole()->partnerDipoles().end(); ++dip )
    sum += channelWeight((**dip).realEmitter(),
			 (**dip).realEmission(),
			 (**dip).realSpectator());
  assert(sum > 0.0);
  return currentChannel / sum;
}

CrossSection MEMatching::dSigHatDR() const {

  assert(realXComb()->lastME2() > 0.0);

  double xme2 = realXComb()->lastME2() * channelWeight();
  xme2 /= 
    pow(dipole()->realEmissionME()->lastXComb().lastAlphaS(),
	dipole()->realEmissionME()->orderInAlphaS());
  xme2 *=
    pow(dipole()->underlyingBornME()->lastXComb().lastAlphaS(),
	dipole()->underlyingBornME()->orderInAlphaS());
  xme2 *= bornPDFWeight(dipole()->underlyingBornME()->lastScale());    

  return
    sqr(hbarc) * 
    realXComb()->jacobian() * 
    subtractionScaleWeight() *
    xme2 /
    (2. * realXComb()->lastSHat());

}

double MEMatching::me2() const {

  assert(bornXComb()->lastME2() > 0.0);

  double bme2 = bornXComb()->lastME2();
  bme2 /=
    pow(dipole()->underlyingBornME()->lastXComb().lastAlphaS(),
	dipole()->underlyingBornME()->orderInAlphaS());

  double rme2 = dipole()->realEmissionME()->me2();
  rme2 /= 
    pow(dipole()->realEmissionME()->lastXComb().lastAlphaS(),
	dipole()->realEmissionME()->orderInAlphaS());
  rme2 *= 
    pow(bornXComb()->lastSHat()/realXComb()->lastSHat(),
	2.*(realCXComb()->mePartonData().size())-8.);

  return
    channelWeight() * (rme2/bme2) *
    splittingScaleWeight();

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MEMatching::persistentOutput(PersistentOStream & ) const {}

void MEMatching::persistentInput(PersistentIStream & , int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MEMatching,Herwig::ShowerApproximation>
  describeHerwigMEMatching("Herwig::MEMatching", "HwMatchbox.so");

void MEMatching::Init() {

  static ClassDocumentation<MEMatching> documentation
    ("MEMatching implements NLO matching with matrix element correction (aka Powheg).");

}

