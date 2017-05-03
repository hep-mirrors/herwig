// -*- C++ -*-
//
// FIqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMsqgxDipole class.
//

#include "FIMsqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FIMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FIMassiveInvertedTildeKinematics.h"

using namespace Herwig;

FIMsqgxDipole::FIMsqgxDipole() 
  : SubtractionDipole() {}

FIMsqgxDipole::~FIMsqgxDipole() {}

IBPtr FIMsqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FIMsqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FIMsqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator < 2 &&
    partons[emission]->id() == ParticleID::g &&
    ((abs(partons[emitter]->id())> 1000000 && abs(partons[emitter]->id())< 1000007) ||
     (abs(partons[emitter]->id())> 2000000 && abs(partons[emitter]->id())< 2000007)) &&
    partons[emitter]->hardProcessMass() != ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FIMsqgxDipole::me2Avg(double ccme2) const {
  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  // extra mass terms cancel
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

//   // NOTE: extra term taken from FIqgxDipole implementation??
//   res *= ( 2./(1.-z+(1.-x)) - 2. +(1.-x)*(1.+3.*x*z) -
//   	   sqr(realEmissionME()->lastXComb().mePartonData()[realEmission()]->hardProcessMass()) / prop * 2.*x);
  // NOTE: CR: extra term switched off in massive implementation for the moment,
  //           mass of realEmission changed to mass of realEmitter
  res *= ( 2./(1.-z+(1.-x)) - 2. -
  	   sqr(realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass()) / prop * 2.*x);

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FIMsqgxDipole::me2() const {
  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  // extra mass terms cancel
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

//   // NOTE: extra term taken from FIqgxDipole implementation
//   res *= ( 2./(1.-z+(1.-x)) -2. +(1.-x)*(1.+3.*x*z) -
//   	   sqr(realEmissionME()->lastXComb().mePartonData()[realEmission()]->hardProcessMass()) / prop * 2.*x);
  // NOTE: CR: extra term switched off in massive implementation for the moment,
  //           mass of realEmission changed to mass of realEmitter
  res *= ( 2./(1.-z+(1.-x)) - 2. -
  	   sqr(realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass()) / prop * 2.*x);

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FIMsqgxDipole::persistentOutput(PersistentOStream &) const {
}

void FIMsqgxDipole::persistentInput(PersistentIStream &, int) {
}

void FIMsqgxDipole::Init() {

  static ClassDocumentation<FIMsqgxDipole> documentation
    ("FIMsqgxDipole");

  DipoleRepository::registerDipole<0,FIMsqgxDipole,FIMassiveTildeKinematics,FIMassiveInvertedTildeKinematics>
    ("FIMsqgxDipole","FIMassiveTildeKinematics","FIMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIMsqgxDipole,SubtractionDipole>
describeHerwigFIMsqgxDipole("Herwig::FIMsqgxDipole", "Herwig.so");
