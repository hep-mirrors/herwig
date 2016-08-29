// -*- C++ -*-
//
// IFggxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFggxDipole class.
//

#include "IFggxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightInvertedTildeKinematics.h"

using namespace Herwig;

IFggxDipole::IFggxDipole() 
  : SubtractionDipole() {}

IFggxDipole::~IFggxDipole() {}

IBPtr IFggxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFggxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFggxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    partons[emitter]->id() == ParticleID::g &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IFggxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res = 1./(1.-x+u) + (1.-x)/x - 1. + x*(1.-x);

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

double IFggxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  
  if ( alpha() < u )
    return 0.0;
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double diag = 1./(1.-x+u)-1.+x*(1.-x);
  Lorentz5Momentum pc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]/u -
    realEmissionME()->lastXComb().meMomenta()[realSpectator()]/(1.-u);

  Energy2 sc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]*
    realEmissionME()->lastXComb().meMomenta()[realSpectator()];
  sc /= u*(1.-u)*(1.-x)/x;

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

void IFggxDipole::persistentOutput(PersistentOStream &) const {
}

void IFggxDipole::persistentInput(PersistentIStream &, int) {
}

void IFggxDipole::Init() {

  static ClassDocumentation<IFggxDipole> documentation
    ("IFggxDipole");

  DipoleRepository::registerDipole<0,IFggxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
    ("IFggxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFggxDipole,SubtractionDipole>
describeHerwigIFggxDipole("Herwig::IFggxDipole", "Herwig.so");
