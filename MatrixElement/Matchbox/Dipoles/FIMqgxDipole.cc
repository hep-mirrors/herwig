// -*- C++ -*-
//
// FIqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMqgxDipole class.
//

#include "FIMqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
//#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightTildeKinematics.h"
//#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightInvertedTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FIMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FIMassiveInvertedTildeKinematics.h"

using namespace Herwig;

FIMqgxDipole::FIMqgxDipole() 
  : SubtractionDipole() {}

FIMqgxDipole::~FIMqgxDipole() {}

IBPtr FIMqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FIMqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FIMqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator < 2 &&
    partons[emission]->id() == ParticleID::g &&
//     abs(partons[emitter]->id()) < 7 &&
    ( abs(partons[emitter]->id()) < 7 || abs(partons[emitter]->id()) == 1000021 ) &&
    partons[emitter]->hardProcessMass() != ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FIMqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  if ( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->id() == 1000021 )
    CF = SM().Nc(); // For the SUSY D_{gluino,g}^a subtraction dipole we need to replace CF by CA=Nc

  // extra mass terms cancel
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  // NOTE: extra term taken from FIqgxDipole implementation
  // NOTE: extra term switched off for the moment in the massive case
  res *= ( 
    2./(1.-z+(1.-x)) - (1.+z) 
    // + (1.-x)*(1.+3.*x*z) 
    - sqr(realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass()) / prop * 2.*x
    );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FIMqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  if ( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->id() == 1000021 )
    CF = SM().Nc(); // For the SUSY D_{gluino,g}^a subtraction dipole we need to replace CF by CA=Nc

  // extra mass terms cancel
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  // NOTE: extra term taken from FIqgxDipole implementation
  // NOTE: extra term switched off for the moment in the massive case
  res *= ( 
    2./(1.-z+(1.-x)) - (1.+z) 
    // + (1.-x)*(1.+3.*x*z) 
    - sqr(realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass()) / prop * 2.*x
    );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FIMqgxDipole::persistentOutput(PersistentOStream &) const {
}

void FIMqgxDipole::persistentInput(PersistentIStream &, int) {
}

void FIMqgxDipole::Init() {

  static ClassDocumentation<FIMqgxDipole> documentation
    ("FIMqgxDipole");

//  DipoleRepository::registerDipole<0,FIMqgxDipole,FILightTildeKinematics,FILightInvertedTildeKinematics>
//    ("FIMqgxDipole","FILightTildeKinematics","FILightInvertedTildeKinematics");
  DipoleRepository::registerDipole<0,FIMqgxDipole,FIMassiveTildeKinematics,FIMassiveInvertedTildeKinematics>
    ("FIMqgxDipole","FIMassiveTildeKinematics","FIMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIMqgxDipole,SubtractionDipole>
describeHerwigFIMqgxDipole("Herwig::FIMqgxDipole", "Herwig.so");
