// -*- C++ -*-
//
// FIMqqxDipole.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIqqxDipole class.
//

#include "FIMqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/FILightTildeKinematics.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/FILightInvertedTildeKinematics.h"

using namespace Herwig;

FIMqqxDipole::FIMqqxDipole() 
  : SubtractionDipole() {}

FIMqqxDipole::~FIMqqxDipole() {}

IBPtr FIMqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FIMqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FIMqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator < 2 &&
    abs(partons[emission]->id()) < 6 &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emission]->id() + partons[emitter]->id() == 0 &&
    !(partons[emitter]->mass() == ZERO &&
      partons[emission]->mass() == ZERO &&
      partons[spectator]->mass() == ZERO);
}

double FIMqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  Energy2 mQ2 = sqr(realEmissionME()->lastXComb().mePartonData()[realEmitter()]->mass());
  double muQ2 = x * mQ2 /
    ((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
     (realEmissionME()->lastXComb().meMomenta()[realSpectator()]));

  // mu_ij=0, mu_i=mu_j=mu_Q.
  double zm = ( 1.-x - sqrt( sqr(1.-x-2.*muQ2) - 4.*muQ2 ) ) /
    ( 2.*(1.-x) );
  double zp = ( 1.-x + sqrt( sqr(1.-x-2.*muQ2) - 4.*muQ2 ) ) /
    ( 2.*(1.-x) );

  double res = 1.-2.*(z-zm)*(zp-z);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/( prop+2.*mQ2*x );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FIMqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  Energy2 mQ2 = sqr((realEmissionME()->lastXComb().mePartonData()[realEmitter()])->mass());

  Lorentz5Momentum pc = 
    z*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    (1.-z)*realEmissionME()->lastXComb().meMomenta()[realEmission()];

  SpinCorrelationTensor corr(-1.,pc,-(prop+2.*mQ2*x)/(4.*x));

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/( prop+2.*mQ2*x );

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FIMqqxDipole::persistentOutput(PersistentOStream &) const {
}

void FIMqqxDipole::persistentInput(PersistentIStream &, int) {
}

void FIMqqxDipole::Init() {

  static ClassDocumentation<FIMqqxDipole> documentation
    ("FIMqqxDipole");

  DipoleRepository::registerDipole<0,FIMqqxDipole,FILightTildeKinematics,FILightInvertedTildeKinematics>
    ("FIMqqxDipole","FILightTildeKinematics","FILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIMqqxDipole,SubtractionDipole>
describeHerwigFIMqqxDipole("Herwig::FIMqqxDipole", "HwMatchbox.so");
