// -*- C++ -*-
//
// FFggxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFggxDipole class.
//

#include "FFggxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFLightInvertedTildeKinematics.h"

using namespace Herwig;

FFggxDipole::FFggxDipole() 
  : SubtractionDipole() {}

FFggxDipole::~FFggxDipole() {}

IBPtr FFggxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFggxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFggxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    partons[emitter]->id() == ParticleID::g &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FFggxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double res = 
    1./(1.-z*(1.-y)) + 1./(1.-(1.-z)*(1.-y)) - 2. + z*(1.-z);

  res *= -ccme2;

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FFggxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  if ( alpha() < y )
    return 0.0;
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double diag = 1./(1.-z*(1.-y))+1./(1.-(1.-z)*(1.-y))-2.;
  Lorentz5Momentum pc = 
    z*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    (1.-z)*realEmissionME()->lastXComb().meMomenta()[realEmission()];

  SpinCorrelationTensor corr(-diag,pc,prop/2.);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FFggxDipole::persistentOutput(PersistentOStream &) const {
}

void FFggxDipole::persistentInput(PersistentIStream &, int) {
}

void FFggxDipole::Init() {

  static ClassDocumentation<FFggxDipole> documentation
    ("FFggxDipole");

  DipoleRepository::registerDipole<0,FFggxDipole,FFLightTildeKinematics,FFLightInvertedTildeKinematics>
    ("FFggxDipole","FFLightTildeKinematics","FFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFggxDipole,SubtractionDipole>
describeHerwigFFggxDipole("Herwig::FFggxDipole", "Herwig.so");
