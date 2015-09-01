// -*- C++ -*-
//
// IFMqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFqgxDipole class.
//

#include "IFMqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

IFMqgxDipole::IFMqgxDipole() 
  : SubtractionDipole() {}

IFMqgxDipole::~IFMqgxDipole() {}

IBPtr IFMqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFMqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFMqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    abs(partons[emitter]->id()) < 7 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() != ZERO;
}

double IFMqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  // NOTE: extra term same as in IFqgxDipole
  // NOTE: extra term switched off for the moment in the massive case
  res *= ( 
    2./(1.-x+u) - (1.+x) 
    // + u*(1.+3.*x*(1.-u)) 
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

double IFMqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  // NOTE: extra term same as in IFqgxDipole
  // NOTE: extra term switched off for the moment in the massive case
  res *= ( 
    2./(1.-x+u) - (1.+x) 
    // + u*(1.+3.*x*(1.-u)) 
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

void IFMqgxDipole::persistentOutput(PersistentOStream &) const {
}

void IFMqgxDipole::persistentInput(PersistentIStream &, int) {
}

void IFMqgxDipole::Init() {

  static ClassDocumentation<IFMqgxDipole> documentation
    ("IFMqgxDipole");

//  DipoleRepository::registerDipole<0,IFMqgxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
//    ("IFMqgxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");
  DipoleRepository::registerDipole<0,IFMqgxDipole,IFMassiveTildeKinematics,IFMassiveInvertedTildeKinematics>
    ("IFMqgxDipole","IFMassiveTildeKinematics","IFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMqgxDipole,SubtractionDipole>
describeHerwigIFMqgxDipole("Herwig::IFMqgxDipole", "Herwig.so");
