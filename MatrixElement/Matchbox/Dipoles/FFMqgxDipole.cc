// -*- C++ -*-
//
// FFMqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMqgxDipole class.
//

#include "FFMqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

FFMqgxDipole::FFMqgxDipole() 
  : SubtractionDipole() {}

FFMqgxDipole::~FFMqgxDipole() {}

IBPtr FFMqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFMqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFMqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
//     abs(partons[emitter]->id()) < 7 &&
    ( abs(partons[emitter]->id()) < 7 || abs(partons[emitter]->id()) == 1000021 ) &&
    !(partons[emitter]->hardProcessMass() == ZERO &&
      partons[spectator]->hardProcessMass() == ZERO);
}

double FFMqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  // masses
  double muQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass() / lastDipoleScale() );
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-muQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-muQ2-muj2)*(1.-y));
  double vbar = sqrt( 1.+sqr(muQ2)+sqr(muj2)-2.*(muQ2+muj2+muQ2*muj2) ) / (1.-muQ2-muj2);

  Energy2 prop =
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
  (realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  if ( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->id() == 1000021 )
    CF = SM().Nc(); // For the SUSY D_{gluino,g;k} subtraction dipole we need to replace CF by CA=Nc

  // extra mass terms cancel: mi2+m2-Mi2 = mQ2+0-mQ2
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z*(1.-y)) - vbar/vijk * ( (1.+z) + muQ2*sqr(lastDipoleScale())*2./prop ) );

  res *= -ccme2;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;
}

double FFMqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  // masses
  double muQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass() / lastDipoleScale() );
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-muQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-muQ2-muj2)*(1.-y));
  double vbar = sqrt( 1.+sqr(muQ2)+sqr(muj2)-2.*(muQ2+muj2+muQ2*muj2) ) / (1.-muQ2-muj2);

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  if ( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->id() == 1000021 )
    CF = SM().Nc(); // For the SUSY D_{gluino,g;k} subtraction dipole we need to replace CF by CA=Nc

  // extra mass terms cancel: mi2+m2-Mi2 = mQ2+0-mQ2
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z*(1.-y)) - vbar/vijk * ( (1.+z) + muQ2*sqr(lastDipoleScale())*2./prop ) );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FFMqgxDipole::persistentOutput(PersistentOStream &) const {
}

void FFMqgxDipole::persistentInput(PersistentIStream &, int) {
}

void FFMqgxDipole::Init() {

  static ClassDocumentation<FFMqgxDipole> documentation
    ("FFMqgxDipole");

  DipoleRepository::registerDipole<0,FFMqgxDipole,FFMassiveTildeKinematics,FFMassiveInvertedTildeKinematics>
    ("FFMqgxDipole","FFMassiveTildeKinematics","FFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMqgxDipole,SubtractionDipole>
describeHerwigFFMqgxDipole("Herwig::FFMqgxDipole", "Herwig.so");
