// -*- C++ -*-
//
// IFqqxDipole.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFqqxDipole class.
//

#include "IFqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/IFLightTildeKinematics.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/IFLightInvertedTildeKinematics.h"

using namespace Herwig;

IFqqxDipole::IFqqxDipole() 
  : SubtractionDipole() {}

IFqqxDipole::~IFqqxDipole() {}

IBPtr IFqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    abs(partons[emission]->id()) < 6 &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emission]->id() - partons[emitter]->id() == 0 &&
    partons[emitter]->mass() == ZERO &&
    partons[emission]->mass() == ZERO &&
    partons[spectator]->mass() == ZERO;
}

double IFqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res = x + 2.*(1.-x)/x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  res *= 8.*CF*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (realEmissionME()->lastXComb().lastAlphaS())/prop;

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

double IFqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  Lorentz5Momentum pc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]/u -
    realEmissionME()->lastXComb().meMomenta()[realSpectator()]/(1.-u);

  Energy2 sc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]*
    realEmissionME()->lastXComb().meMomenta()[realSpectator()];
  sc /= u*(1.-u)*(1.-x)/x;

  SpinCorrelationTensor corr(-x,pc,sc/2.);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  res *= 8.*CF*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (realEmissionME()->lastXComb().lastAlphaS())/prop;

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

void IFqqxDipole::persistentOutput(PersistentOStream &) const {
}

void IFqqxDipole::persistentInput(PersistentIStream &, int) {
}

void IFqqxDipole::Init() {

  static ClassDocumentation<IFqqxDipole> documentation
    ("IFqqxDipole");

  DipoleRepository::registerDipole<IFqqxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
    ("IFqqxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFqqxDipole,SubtractionDipole>
describeHerwigIFqqxDipole("Herwig::IFqqxDipole", "HwMatchbox.so");
