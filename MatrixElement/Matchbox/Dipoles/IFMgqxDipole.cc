// -*- C++ -*-
//
// IFMgqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFgqxDipole class.
//

#include "IFMgqxDipole.h"
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

IFMgqxDipole::IFMgqxDipole() 
  : SubtractionDipole() {}

IFMgqxDipole::~IFMgqxDipole() {}

IBPtr IFMgqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFMgqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFMgqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    partons[emitter]->id() == ParticleID::g &&
    abs(partons[emission]->id()) < 7 &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() != ZERO;
}

double IFMgqxDipole::me2Avg(double ccme2) const {

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

double IFMgqxDipole::me2() const {

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

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IFMgqxDipole::persistentOutput(PersistentOStream &) const {
}

void IFMgqxDipole::persistentInput(PersistentIStream &, int) {
}

void IFMgqxDipole::Init() {

  static ClassDocumentation<IFMgqxDipole> documentation
    ("IFMgqxDipole");

  DipoleRepository::registerDipole<0,IFMgqxDipole,IFMassiveTildeKinematics,IFMassiveInvertedTildeKinematics>
    ("IFMgqxDipole","IFMassiveTildeKinematics","IFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMgqxDipole,SubtractionDipole>
describeHerwigIFMgqxDipole("Herwig::IFMgqxDipole", "Herwig.so");
