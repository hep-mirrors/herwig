// -*- C++ -*-
//
// FIqqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIqqxDipole class.
//

#include "FIqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightInvertedTildeKinematics.h"

using namespace Herwig;

FIqqxDipole::FIqqxDipole() 
  : SubtractionDipole() {}

FIqqxDipole::~FIqqxDipole() {}

IBPtr FIqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FIqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FIqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator < 2 &&
    abs(partons[emission]->id()) < 6 &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emission]->id() + partons[emitter]->id() == 0 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FIqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res = 1.-2.*z*(1.-z);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FIqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  if ( alpha() < (1.-x) )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  Lorentz5Momentum pc = 
    z*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    (1.-z)*realEmissionME()->lastXComb().meMomenta()[realEmission()];

  SpinCorrelationTensor corr(-1.,pc,-prop/(4.*x));

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FIqqxDipole::persistentOutput(PersistentOStream &) const {
}

void FIqqxDipole::persistentInput(PersistentIStream &, int) {
}

void FIqqxDipole::Init() {

  static ClassDocumentation<FIqqxDipole> documentation
    ("FIqqxDipole");

  DipoleRepository::registerDipole<0,FIqqxDipole,FILightTildeKinematics,FILightInvertedTildeKinematics>
    ("FIqqxDipole","FILightTildeKinematics","FILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIqqxDipole,SubtractionDipole>
describeHerwigFIqqxDipole("Herwig::FIqqxDipole", "Herwig.so");
