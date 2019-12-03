// -*- C++ -*-
//
// FFMsqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMsqgxDipole class.
//

#include "FFMsqgxDipole.h"
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

FFMsqgxDipole::FFMsqgxDipole() 
  : SubtractionDipole() {}

FFMsqgxDipole::~FFMsqgxDipole() {}

IBPtr FFMsqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFMsqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFMsqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    ((abs(partons[emitter]->id())> 1000000 && abs(partons[emitter]->id())< 1000007) ||
     (abs(partons[emitter]->id())> 2000000 && abs(partons[emitter]->id())< 2000007)) &&
    !(partons[emitter]->hardProcessMass() == ZERO &&
      partons[spectator]->hardProcessMass() == ZERO);
}


double FFMsqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  // masses
  double muSQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass() / lastDipoleScale() );
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-muSQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-muSQ2-muj2)*(1.-y));
  double vbar = sqrt( 1.+sqr(muSQ2)+sqr(muj2)-2.*(muSQ2+muj2+muSQ2*muj2) ) / (1.-muSQ2-muj2);

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z*(1.-y)) - vbar/vijk * ( 2. + muSQ2*sqr(lastDipoleScale())*2./prop ));

  res *= -ccme2;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;


}

double FFMsqgxDipole::me2() const {
  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  // masses
  double muSQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass() / lastDipoleScale() );
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-muSQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-muSQ2-muj2)*(1.-y));
  double vbar = sqrt( 1.+sqr(muSQ2)+sqr(muj2)-2.*(muSQ2+muj2+muSQ2*muj2) ) / (1.-muSQ2-muj2);

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  // extra mass terms cancel: mi2+m2-Mi2 = mQ2+0-mQ2
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z*(1.-y)) - vbar/vijk * ( 2. + muSQ2*sqr(lastDipoleScale())*2./prop ) );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  
  return res;

}

void FFMsqgxDipole::persistentOutput(PersistentOStream &) const {
}

void FFMsqgxDipole::persistentInput(PersistentIStream &, int) {
}

void FFMsqgxDipole::Init() {

  static ClassDocumentation<FFMsqgxDipole> documentation
    ("FFMsqgxDipole");

  DipoleRepository::registerDipole<0,FFMsqgxDipole,FFMassiveTildeKinematics,FFMassiveInvertedTildeKinematics>
    ("FFMsqgxDipole","FFMassiveTildeKinematics","FFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMsqgxDipole,SubtractionDipole>
describeHerwigFFMsqgxDipole("Herwig::FFMsqgxDipole", "Herwig.so");
