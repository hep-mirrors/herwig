// -*- C++ -*-
//
// IFqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFqgxDipole class.
//

#include "IFqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightInvertedTildeKinematics.h"

using namespace Herwig;

IFqgxDipole::IFqgxDipole() 
  : SubtractionDipole() {}

IFqgxDipole::~IFqgxDipole() {}

IBPtr IFqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IFqgxDipole::me2Avg(double ccme2) const {

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

  res *= ( 
    2./(1.-x+u) - (1.+x) 
    + u*(1.+3.*x*(1.-u)) 
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

double IFqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  
  if ( alpha() < u )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 
    2./(1.-x+u) - (1.+x) 
    //+ u*(1.+3.*x*(1.-u)) 
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

void IFqgxDipole::persistentOutput(PersistentOStream &) const {
}

void IFqgxDipole::persistentInput(PersistentIStream &, int) {
}

void IFqgxDipole::Init() {

  static ClassDocumentation<IFqgxDipole> documentation
    ("IFqgxDipole");

  DipoleRepository::registerDipole<0,IFqgxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
    ("IFqgxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFqgxDipole,SubtractionDipole>
describeHerwigIFqgxDipole("Herwig::IFqgxDipole", "Herwig.so");
