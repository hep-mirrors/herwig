// -*- C++ -*-
//
// IIgqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIgqxDipole class.
//

#include "IIgqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightInvertedTildeKinematics.h"

using namespace Herwig;

IIgqxDipole::IIgqxDipole() 
  : SubtractionDipole() {}

IIgqxDipole::~IIgqxDipole() {}

IBPtr IIgqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IIgqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IIgqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator < 2 &&
    partons[emitter]->id() == ParticleID::g &&
    abs(partons[emission]->id()) < 6 &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IIgqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res =
    8.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= .5 * ( 1.-2.*x*(1.-x) );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IIgqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  
  if ( alpha() < v )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res =
    8.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= .5 * ( 1.-2.*x*(1.-x) );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IIgqxDipole::persistentOutput(PersistentOStream &) const {
}

void IIgqxDipole::persistentInput(PersistentIStream &, int) {
}

void IIgqxDipole::Init() {

  static ClassDocumentation<IIgqxDipole> documentation
    ("IIgqxDipole");

  DipoleRepository::registerDipole<0,IIgqxDipole,IILightTildeKinematics,IILightInvertedTildeKinematics>
    ("IIgqxDipole","IILightTildeKinematics","IILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IIgqxDipole,SubtractionDipole>
describeHerwigIIgqxDipole("Herwig::IIgqxDipole", "Herwig.so");
