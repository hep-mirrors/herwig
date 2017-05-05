// -*- C++ -*-
//
// FFMassiveTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMassiveTildeKinematics class.
//

#include "FFMassiveTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFMassiveTildeKinematics::FFMassiveTildeKinematics() {}

FFMassiveTildeKinematics::~FFMassiveTildeKinematics() {}

IBPtr FFMassiveTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFMassiveTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FFMassiveTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  // Compute y
  double y = emission*emitter / (emission*emitter + emission*spectator + emitter*spectator);

  // Calculate the scale
  Lorentz5Momentum pTot = emitter+emission+spectator;
  Energy scale = pTot.m();

  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  double Mui2 = sqr( bornEmitterData()->hardProcessMass() / scale );
  double Muj2 = sqr( bornSpectatorData()->hardProcessMass() / scale );
 
  // Calculate the invariants
  Energy2 Qijk = sqr(scale);
  double suijk = 0.5*( 1. - Mui2 - Muj2 + sqrt( sqr(1.-Mui2-Muj2) - 4.*Mui2*Muj2 ) );
  double bar = 1. - mui2 - mu2 - muj2;

  // Calculate A (as in notes)
  double A = (y*bar - Mui2 + mui2 + mu2)/suijk;

  // Calculate the scale factors, xk and xij
  double lambdaK = 1. + (Muj2/suijk);
  double lambdaIJ = 1. + (Mui2/suijk);
  double xk = (1./(2.*lambdaK)) * ( (lambdaK + (Muj2/suijk)*lambdaIJ - A) + sqrt( sqr(lambdaK + (Muj2/suijk)*lambdaIJ - A) - 4.*lambdaK*lambdaIJ*Muj2/suijk) );
  double xij = 1. - ( (Muj2/suijk) * (1.-xk) / xk );

  // Calculate zPrime = qi.nk / (qi+qj).nk from qi.qk and y
  double l = Muj2 / (2.*xk*xij*suijk);
  double a = (xij*xk*suijk/2.) - l*(bar*y + mui2 + mu2);
  double b = (-emitter*spectator)/Qijk + l*(bar*y + 2.*mui2);
  double zPrime = -b/a;


  // Store z as well
  double z = emitter*spectator / ((emitter + emission)*spectator);

  subtractionParameters().resize(3);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = z;
  subtractionParameters()[2] = zPrime;
  
  // Calculate nij and nk from qi, q, qj
  double B = (Mui2/(xij*suijk) + A/xij) / xk;
  double den = (xij/B - Muj2/(xk*suijk));

  Lorentz5Momentum nij = (1./den)*( (emitter+emission)/B - spectator);
  Lorentz5Momentum nk = (1./xk)*(spectator - Muj2*nij/(xk*suijk));

  // Calculate the born momenta from nij and nk
  bornSpectatorMomentum() = nk + (Muj2/suijk)*nij;
  bornEmitterMomentum() = pTot - bornSpectatorMomentum();

  bornEmitterMomentum().setMass( bornEmitterData()->hardProcessMass() );
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass( bornSpectatorData()->hardProcessMass() );
  bornSpectatorMomentum().rescaleEnergy();
  
  return true;

}

Energy FFMassiveTildeKinematics::lastPt() const {
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
  
  double y = subtractionParameters()[0];
  double zPrime = subtractionParameters()[2];
  
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / scale );

  Energy ret = scale * sqrt( y * (1.-mui2-mu2-muj2) * zPrime*(1.-zPrime) - sqr(1.-zPrime)*mui2 - sqr(zPrime)*mu2 );

  return ret;
}


Energy FFMassiveTildeKinematics::lastPt(Lorentz5Momentum emitter,
					Lorentz5Momentum emission,
					Lorentz5Momentum spectator)const {


  // Compute y
  double y = emission*emitter / (emission*emitter 
			       + emission*spectator 
			       + emitter*spectator);

  // Calculate the scale
  Lorentz5Momentum pTot = emitter+emission+spectator;
  Energy scale = pTot.m();

  // masses
  double mui2 = sqr( emitter.mass() / scale );
  double mu2  = sqr( emission.mass() / scale );
  double muj2 = sqr( spectator.mass() / scale );
   // TODO: here we assume a gluon
  bool isgluon= emitter.mass()==emission.mass();
  double Mui2 = sqr(( isgluon?ZERO:max(emission.mass(),emitter.mass()) )/ scale );
  double Muj2 = sqr( spectator.mass() / scale );
 
  // Calculate the invariants
  Energy2 Qijk = sqr(scale);
  double suijk = 0.5*( 1. - Mui2 - Muj2 + sqrt( sqr(1.-Mui2-Muj2) - 4.*Mui2*Muj2 ) );
  double bar = 1. - mui2 - mu2 - muj2;

  // Calculate A (as in notes)
  double A = (y*bar - Mui2 + mui2 + mu2)/suijk;

  // Calculate the scale factors, xk and xij
  double lambdaK = 1. + (Muj2/suijk);
  double lambdaIJ = 1. + (Mui2/suijk);
  double xk = (1./(2.*lambdaK)) * ( (lambdaK + (Muj2/suijk)*lambdaIJ - A) + sqrt( sqr(lambdaK + (Muj2/suijk)*lambdaIJ - A) - 4.*lambdaK*lambdaIJ*Muj2/suijk) );
  double xij = 1. - ( (Muj2/suijk) * (1.-xk) / xk );

  // Calculate zPrime = qi.nk / (qi+qj).nk from qi.qk and y
  double l = Muj2 / (2.*xk*xij*suijk);
  double a = (xij*xk*suijk/2.) - l*(bar*y + mui2 + mu2);
  double b = (-emitter*spectator)/Qijk + l*(bar*y + 2.*mui2);
  double zPrime = -b/a;



  Energy ret = scale * sqrt( y * (1.-mui2-mu2-muj2) * zPrime*(1.-zPrime) - sqr(1.-zPrime)*mui2 - sqr(zPrime)*mu2 );

  return ret;


}


  // NOTE: bounds calculated at this step may be too loose
pair<double,double> FFMassiveTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  
  if(pt>hardPt) return make_pair(0.5,0.5);
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
    // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  
  double zp = ( 1.+mui2-mu2+muj2-2.*sqrt(muj2) +
               rootOfKallen(mui2,mu2,sqr(1-sqrt(muj2))) *
               sqrt( 1.-sqr(pt/hardPt) ) ) /
  ( 2.*sqr(1.-sqrt(muj2)) );
  double zm = ( 1.+mui2-mu2+muj2-2.*sqrt(muj2) -
               rootOfKallen(mui2,mu2,sqr(1-sqrt(muj2))) *
               sqrt( 1.-sqr(pt/hardPt) ) ) /
  ( 2.*sqr(1.-sqrt(muj2)) );
  
  return make_pair(zm,zp);
}



double FFMassiveTildeKinematics::lastZ() const {
  return subtractionParameters()[2];
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMassiveTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FFMassiveTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FFMassiveTildeKinematics::Init() {

  static ClassDocumentation<FFMassiveTildeKinematics> documentation
    ("FFMassiveTildeKinematics implements the 'tilde' kinematics for "
     "a final-final subtraction dipole involving a massive particle.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMassiveTildeKinematics,TildeKinematics>
describeHerwigFFMassiveTildeKinematics("Herwig::FFMassiveTildeKinematics", "Herwig.so");
