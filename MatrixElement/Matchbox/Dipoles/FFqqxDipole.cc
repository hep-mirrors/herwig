// -*- C++ -*-
//
// FFqqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFqqxDipole class.
//

#include "FFqqxDipole.h"
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

FFqqxDipole::FFqqxDipole() 
  : SubtractionDipole() {}

FFqqxDipole::~FFqqxDipole() {}

IBPtr FFqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    abs(partons[emission]->id()) < 6 &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emission]->id() + partons[emitter]->id() == 0 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FFqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double z = subtractionParameters()[1];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double res = 1.-2.*z*(1.-z);

  res *= -ccme2;

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FFqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;
  
  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
 
  if ( alpha() < y )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  Lorentz5Momentum pc = 
    z*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    (1.-z)*realEmissionME()->lastXComb().meMomenta()[realEmission()];

  SpinCorrelationTensor corr(-1.,pc,-prop/4.);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FFqqxDipole::persistentOutput(PersistentOStream &) const {
}

void FFqqxDipole::persistentInput(PersistentIStream &, int) {
}

void FFqqxDipole::Init() {

  static ClassDocumentation<FFqqxDipole> documentation
    ("FFqqxDipole");

  DipoleRepository::registerDipole<0,FFqqxDipole,FFLightTildeKinematics,FFLightInvertedTildeKinematics>
    ("FFqqxDipole","FFLightTildeKinematics","FFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFqqxDipole,SubtractionDipole>
describeHerwigFFqqxDipole("Herwig::FFqqxDipole", "Herwig.so");
