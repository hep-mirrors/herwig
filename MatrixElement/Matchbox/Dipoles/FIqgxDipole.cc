// -*- C++ -*-
//
// FIqgxDipole.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIqgxDipole class.
//

#include "FIqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/FILightTildeKinematics.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/FILightInvertedTildeKinematics.h"

using namespace Herwig;

FIqgxDipole::FIqgxDipole() 
  : SubtractionDipole() {}

FIqgxDipole::~FIqgxDipole() {}

IBPtr FIqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FIqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FIqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator < 2 &&
    partons[emission]->id() == ParticleID::g &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emitter]->mass() == ZERO &&
    partons[spectator]->mass() == ZERO;
}

double FIqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (realEmissionME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z+(1.-x)) -(1.+z) +(1.-x)*(1.+3.*x*z) );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

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

  return res;

}

double FIqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (realEmissionME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z+(1.-x)) -(1.+z) +(1.-x)*(1.+3.*x*z) );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

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

void FIqgxDipole::persistentOutput(PersistentOStream &) const {
}

void FIqgxDipole::persistentInput(PersistentIStream &, int) {
}

void FIqgxDipole::Init() {

  static ClassDocumentation<FIqgxDipole> documentation
    ("FIqgxDipole");

  DipoleRepository::registerDipole<FIqgxDipole,FILightTildeKinematics,FILightInvertedTildeKinematics>
    ("FIqgxDipole","FILightTildeKinematics","FILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIqgxDipole,SubtractionDipole>
describeHerwigFIqgxDipole("Herwig::FIqgxDipole", "HwMatchbox.so");
