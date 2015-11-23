// -*- C++ -*-
//
// IFMqqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFqqxDipole class.
//

#include "IFMqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
//#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightTildeKinematics.h"
//#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightInvertedTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

IFMqqxDipole::IFMqqxDipole() 
  : SubtractionDipole() {}

IFMqqxDipole::~IFMqqxDipole() {}

IBPtr IFMqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFMqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFMqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
//    abs(partons[emission]->id()) < 6 &&
//    abs(partons[emitter]->id()) < 6 &&
    abs(partons[emission]->id()) < 7 &&
    abs(partons[emitter]->id()) < 7 &&
    partons[emission]->id() - partons[emitter]->id() == 0 &&
//    !(partons[emitter]->hardProcessMass() == ZERO &&
//      partons[emission]->hardProcessMass() == ZERO &&
//      partons[spectator]->hardProcessMass() == ZERO);
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() != ZERO;
}

double IFMqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double muj2 = sqr( (realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass()) ) /
    (2.* (realEmissionME()->lastXComb().meMomenta()[bornSpectator()])*
     (realEmissionME()->lastXComb().meMomenta()[realEmitter()]) );

  double res = x + 2.*(1.-x)/x -
    2.*muj2/x*u/(1.-u);

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  res *= 8.*CF*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
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

double IFMqqxDipole::me2() const {

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
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IFMqqxDipole::persistentOutput(PersistentOStream &) const {
}

void IFMqqxDipole::persistentInput(PersistentIStream &, int) {
}

void IFMqqxDipole::Init() {

  static ClassDocumentation<IFMqqxDipole> documentation
    ("IFMqqxDipole");

//  DipoleRepository::registerDipole<0,IFMqqxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
//    ("IFMqqxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");
  DipoleRepository::registerDipole<0,IFMqqxDipole,IFMassiveTildeKinematics,IFMassiveInvertedTildeKinematics>
    ("IFMqqxDipole","IFMassiveTildeKinematics","IFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMqqxDipole,SubtractionDipole>
describeHerwigIFMqqxDipole("Herwig::IFMqqxDipole", "Herwig.so");
