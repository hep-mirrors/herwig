// -*- C++ -*-
//
// FFMggxDipole.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMggxDipole class.
//

#include "FFMggxDipole.h"
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

FFMggxDipole::FFMggxDipole() 
  : SubtractionDipole() {}

FFMggxDipole::~FFMggxDipole() {}

IBPtr FFMggxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFMggxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFMggxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    partons[emitter]->id() == ParticleID::g &&
    partons[spectator]->mass() != ZERO;
}

// TODO: (5.20)
double FFMggxDipole::me2Avg(double) const {
  assert(false && "implementation missing");
  return 0.;
}

double FFMggxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  // masses, g->gg all masses zero except spectator
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->mass() / lastDipoleScale() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-muj2)*(1.-y))-4.*muj2 ) / ((1.-muj2)*(1.-y));
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double diag = 1./(1-z*(1-y))+1./(1-(1.-z)*(1.-y))-2./vijk; // kappa=0
  double zim = z-0.5*(1.-vijk), zjm = (1.-z)-0.5*(1.-vijk);
  Lorentz5Momentum pc = 
    zim*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    zjm*realEmissionME()->lastXComb().meMomenta()[realEmission()];
  
  SpinCorrelationTensor corr(-diag,pc,prop/2.*vijk);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  // extra mass terms all = 0.
  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (realEmissionME()->lastXComb().lastAlphaS())/prop;

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

void FFMggxDipole::persistentOutput(PersistentOStream &) const {
}

void FFMggxDipole::persistentInput(PersistentIStream &, int) {
}

void FFMggxDipole::Init() {

  static ClassDocumentation<FFMggxDipole> documentation
    ("FFMggxDipole");

  DipoleRepository::registerDipole<FFMggxDipole,FFMassiveTildeKinematics,FFMassiveInvertedTildeKinematics>
    ("FFMggxDipole","FFMassiveTildeKinematics","FFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMggxDipole,SubtractionDipole>
describeHerwigFFMggxDipole("Herwig::FFMggxDipole", "HwMatchbox.so");
