// -*- C++ -*-
//
// IIqqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIqqxDipole class.
//

#include "IIqqxDipole.h"
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

IIqqxDipole::IIqqxDipole() 
  : SubtractionDipole() {}

IIqqxDipole::~IIqqxDipole() {}

IBPtr IIqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IIqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IIqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator < 2 &&
    abs(partons[emission]->id()) < 6 &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emission]->id() - partons[emitter]->id() == 0 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IIqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res = (1.+sqr(1.-x))/x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  res *= 8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
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

double IIqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  
  if ( alpha() < v )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  Lorentz5Momentum pc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()] -
    v*realEmissionME()->lastXComb().meMomenta()[realSpectator()];

  Energy2 sc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]*
    realEmissionME()->lastXComb().meMomenta()[realSpectator()];
  sc /= (1.-x)/(x*v);

  SpinCorrelationTensor corr(-x,pc,sc/2.);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  res *= 8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IIqqxDipole::persistentOutput(PersistentOStream &) const {
}

void IIqqxDipole::persistentInput(PersistentIStream &, int) {
}

void IIqqxDipole::Init() {

  static ClassDocumentation<IIqqxDipole> documentation
    ("IIqqxDipole");

  DipoleRepository::registerDipole<0,IIqqxDipole,IILightTildeKinematics,IILightInvertedTildeKinematics>
    ("IIqqxDipole","IILightTildeKinematics","IILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IIqqxDipole,SubtractionDipole>
describeHerwigIIqqxDipole("Herwig::IIqqxDipole", "Herwig.so");
