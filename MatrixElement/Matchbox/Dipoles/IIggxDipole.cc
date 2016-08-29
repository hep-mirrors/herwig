// -*- C++ -*-
//
// IIggxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIggxDipole class.
//

#include "IIggxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightInvertedTildeKinematics.h"

using namespace Herwig;

IIggxDipole::IIggxDipole() 
  : SubtractionDipole() {}

IIggxDipole::~IIggxDipole() {}

IBPtr IIggxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IIggxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IIggxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator < 2 &&
    partons[emission]->id() == ParticleID::g &&
    partons[emitter]->id() == ParticleID::g &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IIggxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res = x/(1.-x) + (1.-x)/x + x*(1.-x);

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
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

double IIggxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  
  if ( alpha() < v )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double diag = x/(1.-x)+x*(1.-x);
  Lorentz5Momentum pc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()] -
    v*realEmissionME()->lastXComb().meMomenta()[realSpectator()];

  Energy2 sc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]*
    realEmissionME()->lastXComb().meMomenta()[realSpectator()];
  sc /= (1.-x)/(x*v);

  SpinCorrelationTensor corr(-diag,pc,sc);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IIggxDipole::persistentOutput(PersistentOStream &) const {
}

void IIggxDipole::persistentInput(PersistentIStream &, int) {
}

void IIggxDipole::Init() {

  static ClassDocumentation<IIggxDipole> documentation
    ("IIggxDipole");

  DipoleRepository::registerDipole<0,IIggxDipole,IILightTildeKinematics,IILightInvertedTildeKinematics>
    ("IIggxDipole","IILightTildeKinematics","IILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IIggxDipole,SubtractionDipole>
describeHerwigIIggxDipole("Herwig::IIggxDipole", "Herwig.so");
