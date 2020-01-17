// -*- C++ -*-
//
// MEMatching.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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

#include "ThePEG/Handlers/StdXCombGroup.h"
#include "ThePEG/PDT/EnumParticles.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

using namespace Herwig;

MEMatching::MEMatching()
  : theTruncatedShower(false) {}

MEMatching::~MEMatching() {}

IBPtr MEMatching::clone() const {
  return new_ptr(*this);
}

IBPtr MEMatching::fullclone() const {
  return new_ptr(*this);
}

CrossSection MEMatching::dSigHatDR() const {

  // need to ensure the partner dipoles all got their xcombs set
  // for a safe evaluation of channelweight
  Ptr<StdXCombGroup>::tcptr grp =
    dynamic_ptr_cast<Ptr<StdXCombGroup>::tcptr>(realCXComb());
  assert(grp);
  for ( vector<StdXCombPtr>::const_iterator dep = grp->dependent().begin();
	dep != grp->dependent().end(); ++dep ) {
    (**dep).matrixElement()->setXComb(*dep);
  }

  double xme2 = dipole()->realEmissionME()->me2();

  xme2 *= channelWeight();
  xme2 /= 
    pow(dipole()->realEmissionME()->lastXComb().lastAlphaS(),
	(double)(dipole()->realEmissionME()->orderInAlphaS()));
  xme2 *=
    pow(dipole()->underlyingBornME()->lastXComb().lastAlphaS(),
	(double)(dipole()->underlyingBornME()->orderInAlphaS()));
  double bornPDF = bornPDFWeight(dipole()->underlyingBornME()->lastScale());
  if ( bornPDF == 0.0 )
    return ZERO;

  xme2 *= bornPDF;

  if ( profileScales() )
    xme2 *= profileScales()->hardScaleProfile(dipole()->showerHardScale(),dipole()->lastPt());

  return
    sqr(hbarc) * 
    realXComb()->jacobian() * 
    subtractionScaleWeight() *
    xme2 /
    (2. * realXComb()->lastSHat());

}

double MEMatching::me2() const {

  double bme2 = bornXComb()->matrixElement()->me2();
  bme2 /=
    pow(dipole()->underlyingBornME()->lastXComb().lastAlphaS(),
	(double)(dipole()->underlyingBornME()->orderInAlphaS()));

  double rme2 = dipole()->realEmissionME()->me2();
  rme2 /= 
    pow(dipole()->realEmissionME()->lastXComb().lastAlphaS(),
	(double)(dipole()->realEmissionME()->orderInAlphaS()));
  rme2 *= 
    pow(bornXComb()->lastSHat()/realXComb()->lastSHat(),
	realCXComb()->mePartonData().size()-4.);

  if ( profileScales() )
    rme2 *= profileScales()->hardScaleProfile(dipole()->showerHardScale(),dipole()->lastPt());

  return
    channelWeight() * (rme2/bme2) * 
    (bornXComb()->lastSHat()/realXComb()->lastSHat()) *
    splittingScaleWeight();

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MEMatching::persistentOutput(PersistentOStream & os) const {
  os << theTruncatedShower;
}

void MEMatching::persistentInput(PersistentIStream & is, int) {
  is >> theTruncatedShower;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MEMatching,Herwig::ShowerApproximation>
  describeHerwigMEMatching("Herwig::MEMatching", "Herwig.so");

void MEMatching::Init() {

  static ClassDocumentation<MEMatching> documentation
    ("MEMatching implements NLO matching with matrix element correction (aka Powheg).");

  static Switch<MEMatching,bool> interfaceTruncatedShower
    ("TruncatedShower",
     "Include a truncated qtilde shower",
     &MEMatching::theTruncatedShower, false, false, false);
  static SwitchOption interfaceTruncatedShowerYes
    (interfaceTruncatedShower,
     "Yes",
     "",
     true);
  static SwitchOption interfaceTruncatedShowerNo
    (interfaceTruncatedShower,
     "No",
     "",
     false);

}

