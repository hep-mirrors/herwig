// -*- C++ -*-
//
// FFMqqxDipole.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMqqxDipole class.
//

#include "FFMqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/FFMassiveTildeKinematics.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/FFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

FFMqqxDipole::FFMqqxDipole() 
  : SubtractionDipole() {}

FFMqqxDipole::~FFMqqxDipole() {}

IBPtr FFMqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFMqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFMqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    abs(partons[emission]->id()) < 6 &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emission]->id() + partons[emitter]->id() == 0 &&
    !(partons[emission]->mass() == ZERO &&
      partons[emitter]->mass() == ZERO &&
      partons[spectator]->mass() == ZERO);
}

// TODO: (5.18)
double FFMqqxDipole::me2Avg(double) const {
  assert(false && "implementation missing");
  return 0.;
}

double FFMqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  // masses
  double muQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmission()]->mass() / lastDipoleScale() );
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->mass() / lastDipoleScale() );
  Energy2 mQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmission()]->mass() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-2.*muQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-2.*muQ2-muj2)*(1.-y));

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double zim = z-0.5*(1.-vijk), zjm = (1.-z)-0.5*(1.-vijk);
  Lorentz5Momentum pc = 
    zim*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    zjm*realEmissionME()->lastXComb().meMomenta()[realEmission()];

  // kappa=0 -- otherwise: extra diagonal term (instead of just -1.)
  SpinCorrelationTensor corr(-1.,pc,-(prop+2.*mQ2)/4.);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (realEmissionME()->lastXComb().lastAlphaS())/ (prop+2.*mQ2);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  res *=
    pow(realEmissionME()->lastXComb().lastAlphaS()/
	underlyingBornME()->lastXComb().lastAlphaS(),
	underlyingBornME()->orderInAlphaS());

  res *=
    pow(realEmissionME()->lastXComb().lastAlphaEM()/
	underlyingBornME()->lastXComb().lastAlphaEM(),
	underlyingBornME()->orderInAlphaEW());

  lastME2(res);

  logME2();
  
  return res;

}

void FFMqqxDipole::persistentOutput(PersistentOStream &) const {
}

void FFMqqxDipole::persistentInput(PersistentIStream &, int) {
}

void FFMqqxDipole::Init() {

  static ClassDocumentation<FFMqqxDipole> documentation
    ("FFMqqxDipole");

  DipoleRepository::registerDipole<FFMqqxDipole,FFMassiveTildeKinematics,FFMassiveInvertedTildeKinematics>
    ("FFMqqxDipole","FFMassiveTildeKinematics","FFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMqqxDipole,SubtractionDipole>
describeHerwigFFMqqxDipole("Herwig::FFMqqxDipole", "HwMatchbox.so");
