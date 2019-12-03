// -*- C++ -*-
//
// FFMassiveTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
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

// Matches Stephen Webster's thesis
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
  double muj2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  double Muij2 = sqr( bornEmitterData()->hardProcessMass() / scale );
  double Muk2 = sqr( bornSpectatorData()->hardProcessMass() / scale );
 
  // Calculate the invariants
  Energy2 Qijk = sqr(scale);
  double sijkN = 0.5*( 1. - Muij2 - Muk2 + sqrt( sqr(1.-Muij2-Muk2) - 4.*Muij2*Muk2 ) );
  double bar = 1. - mui2 - muj2 - muk2;
  
  // Calculate Qij2
  double QijN2 = sqr(emitter + emission)/Qijk;
  
  // Calculate the scale factors, xk and xij
  double lambdaK = 1. + (Muk2/sijkN);
  double lambdaIJ = 1. + (Muij2/sijkN);
  double fac1 = lambdaIJ*lambdaK + (muk2 - QijN2)/sijkN;
  double xk =
    ( fac1 + sqrt( sqr(fac1) - 4.*lambdaIJ*lambdaK*muk2/sijkN ) )
    / 2. / lambdaK ;
  double xij = 1. - muk2*(1.-xk) / xk / sijkN;
  
  // Calculate z = qi.nk / (qi+qj).nk from qi.qk and y
  double l = Muk2 / (2.*xk*xij*sijkN);
  double a = (xij*xk*sijkN/2.) - l*(bar*y + mui2 + muj2);
  double b = (-emitter*spectator)/Qijk + l*(bar*y + 2.*mui2);
  double z = -b/a;

  // Calculate zi
  double zi = emitter*spectator / ((emitter + emission)*spectator);

  // Store the variables
  subtractionParameters().resize(3);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = zi;
  subtractionParameters()[2] = z;
  
  // Calculate nij and nk from qi, q, qj
  double L = QijN2/xij/xk/sijkN;
  Lorentz5Momentum nij = (emitter + emission - L*spectator)/(xij - L*Muk2/xk/sijkN);
  Lorentz5Momentum nk = (1./xk)*(spectator - Muk2*nij/(xk*sijkN));
    
  // Calculate the born momenta from nij and nk
  bornSpectatorMomentum() = nk + (Muk2/sijkN)*nij;
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
  double z = subtractionParameters()[2];
  
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / scale );

  Energy ret = scale * sqrt( y * (1.-mui2-muj2-muk2) * z*(1.-z) - sqr(1.-z)*mui2 - sqr(z)*muj2 );
  
  return ret;
  
}

// Matches Stephen Webster's thesis
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
  double muj2  = sqr( emission.mass() / scale );
  double muk2 = sqr( spectator.mass() / scale );
   // TODO: here we assume a gluon
  bool isgluon= emitter.mass()==emission.mass();
  double Muij2 = sqr(( isgluon?ZERO:max(emission.mass(),emitter.mass()) )/ scale );
  double Muk2 = sqr( spectator.mass() / scale );
 
  // Calculate the invariants
  Energy2 Qijk = sqr(scale);
  double sijkN = 0.5*( 1. - Muij2 - Muk2 + sqrt( sqr(1.-Muij2-Muk2) - 4.*Muij2*Muk2 ) );
  double bar = 1. - mui2 - muj2 - muk2;
  
  // Calculate Qij2
  double QijN2 = sqr(emitter + emission)/Qijk;
  
  // Calculate the scale factors, xk and xij
  double lambdaK = 1. + (Muk2/sijkN);
  double lambdaIJ = 1. + (Muij2/sijkN);
  double fac1 = lambdaIJ*lambdaK + (muk2 - QijN2)/sijkN;
  double xk =
    ( fac1 + sqrt( sqr(fac1) - 4.*lambdaIJ*lambdaK*muk2/sijkN ) )
    / 2. / lambdaK ;
  double xij = 1. - muk2*(1.-xk) / xk / sijkN;
  
  // Calculate z = qi.nk / (qi+qj).nk from qi.qk and y
  double l = Muk2 / (2.*xk*xij*sijkN);
  double a = (xij*xk*sijkN/2.) - l*(bar*y + mui2 + muj2);
  double b = (-emitter*spectator)/Qijk + l*(bar*y + 2.*mui2);
  double z = -b/a;

  Energy ret = scale * sqrt( z*(1.-z)*QijN2 - (1.-z)*mui2 - z*muj2 );
    
  return ret;
  
}


  // NOTE: bounds calculated at this step may be too loose
pair<double,double> FFMassiveTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  
  if(pt>hardPt) return make_pair(0.5,0.5);
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
    // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  
  double zp = ( 1.+mui2-muj2+muk2-2.*sqrt(muk2) +
               rootOfKallen(mui2,muj2,sqr(1.-sqrt(muk2))) *
               sqrt( 1.-sqr(pt/hardPt) ) ) /
  ( 2.*sqr(1.-sqrt(muk2)) );
  double zm = ( 1.+mui2-muj2+muk2-2.*sqrt(muk2) -
               rootOfKallen(mui2,muj2,sqr(1.-sqrt(muk2))) *
               sqrt( 1.-sqr(pt/hardPt) ) ) /
  ( 2.*sqr(1.-sqrt(muk2)) );
  
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
