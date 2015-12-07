// -*- C++ -*-
//
// FFMassiveInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMassiveInvertedTildeKinematics class.
//

#include "FFMassiveInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

using namespace Herwig;

FFMassiveInvertedTildeKinematics::FFMassiveInvertedTildeKinematics() {}

FFMassiveInvertedTildeKinematics::~FFMassiveInvertedTildeKinematics() {}

IBPtr FFMassiveInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFMassiveInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FFMassiveInvertedTildeKinematics::doMap(const double * r) {

  if ( ptMax() < ptCut() ) {
    jacobian(0.0);
    return false;
  }

  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();

  double mapping = 1.0;
  pair<Energy,double> ptz = generatePtZ(mapping,r);
  if ( mapping == 0.0 ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;
  double z = ptz.second;
  
  Energy scale = (emitter+spectator).m();
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  double Mui2 = sqr( bornEmitterData()->hardProcessMass() / scale );
  double Muj2 = sqr( bornSpectatorData()->hardProcessMass() / scale );
  
  double y = ( sqr( pt / scale ) + sqr(1.-z)*mui2 + z*z*mu2 ) /
      (z*(1.-z)*(1.-mui2-mu2-muj2));
  
  // now this is done in ptzAllowed!!
  // check (y,z) phasespace boundary
  // 2011-11-09: never occurred
  /**
  double bar = 1.-mui2-mu2-muj2;
  double ym = 2.*sqrt(mui2)*sqrt(mu2)/bar;
  double yp = 1. - 2.*sqrt(muj2)*(1.-sqrt(muj2))/bar;
  double zm = ( (2.*mui2+bar*y)*(1.-y) - sqrt(y*y-ym*ym)*sqrt(sqr(2.*muj2+bar-bar*y)-4.*muj2) ) /
    ( 2.*(1.-y)*(mui2+mu2+bar*y) );
  double zp = ( (2.*mui2+bar*y)*(1.-y) + sqrt(y*y-ym*ym)*sqrt(sqr(2.*muj2+bar-bar*y)-4.*muj2) ) /
    ( 2.*(1.-y)*(mui2+mu2+bar*y) );
  if ( y<ym || y>yp || z<zm || z>zp ) {
    cout << "A problem occurred in FFMassiveInvertedTildeKinematics::doMap: ";
    if(y<ym) cout << "y<ym " << endl;
    if(y>yp) cout << "y>yp " << endl;
    if(z<zm) cout << "z<zm " << endl;
    if(z>zp) cout << "z>zp " << endl;
    jacobian(0.0);
    return false;
  }
  **/
      
  mapping /= z*(1.-z);
  jacobian( mapping*(1.-y)*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)) *
    sqr(1.-mui2-mu2-muj2) / rootOfKallen(1.,Mui2,Muj2) );

  double phi = 2.*Constants::pi*r[2];
  /* // not used ???
  Lorentz5Momentum kt
    = getKt(emitter,spectator,pt,phi);
  */

  subtractionParameters().resize(2);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = z;
  
  // kinematics from FFMassiveKinematics.cc
  // updated 2011-08-23
  // updated 2011-11-08
  Energy2 sbar = sqr(scale) * (1.-mui2-mu2-muj2);
  // CMF: particle energies
  Energy Ei = ( sbar*(1.-(1.-z)*(1.-y)) + 2.*sqr(scale)*mui2 ) / (2.*scale);
  Energy E  = ( sbar*(1.-    z *(1.-y)) + 2.*sqr(scale)*mu2  ) / (2.*scale);
  Energy Ej = ( sbar*(1.-           y ) + 2.*sqr(scale)*muj2 ) / (2.*scale);
  // CMF: momenta in z-direction (axis of pEmitter &pEmitter pSpectator)  
  Energy qi3 = (2.*Ei*Ej-z*(1.-y)*sbar     ) / 2./sqrt(Ej*Ej-sqr(scale)*muj2);
  Energy q3  = (2.*E *Ej-(1.-z)*(1.-y)*sbar) / 2./sqrt(Ej*Ej-sqr(scale)*muj2);
  Energy qj3 = sqrt(   sqr(Ej) - sqr(scale)*muj2 );
  // get z axis in the dipole's CMF which is parallel to pSpectator
  Boost toCMF = (emitter+spectator).findBoostToCM();
  Lorentz5Momentum pjAux = spectator; pjAux.boost(toCMF);
  ThreeVector<double> pjAxis = pjAux.vect().unit();
  // set the momenta in this special reference frame
  // note that pt might in some cases differ from the physical pt!
  Energy ptResc = sqrt( sqr(Ei)-sqr(scale)*mui2-sqr(qi3) );
  Lorentz5Momentum em  ( ptResc*cos(phi), -ptResc*sin(phi), qi3, Ei );
  Lorentz5Momentum emm ( -ptResc*cos(phi), ptResc*sin(phi), q3, E );
  Lorentz5Momentum spe ( 0.*GeV, 0.*GeV, qj3, Ej );
  // rotate back
  em.rotateUz (pjAxis);
  emm.rotateUz(pjAxis);
  spe.rotateUz(pjAxis);
  // boost back
  em.boost (-toCMF);
  emm.boost(-toCMF);
  spe.boost(-toCMF);
  // mass shells, rescale energy
  em.setMass(scale*sqrt(mui2));
  em.rescaleEnergy();
  emm.setMass(scale*sqrt(mu2));
  emm.rescaleEnergy();
  spe.setMass(scale*sqrt(muj2));
  spe.rescaleEnergy();

  // book
  realEmitterMomentum() = em;
  realEmissionMomentum() = emm;
  realSpectatorMomentum() = spe;
  
  return true;

}

Energy FFMassiveInvertedTildeKinematics::lastPt() const {
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  
  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  Energy ret = scale * sqrt( y * (1.-mui2-mu2-muj2) * z*(1.-z) - sqr(1.-z)*mui2 - sqr(z)*mu2 );

  return ret;
}

double FFMassiveInvertedTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}
    
Energy FFMassiveInvertedTildeKinematics::ptMax() const {
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  
  Energy ptmax = rootOfKallen( mui2, mu2, sqr(1.-sqrt(muj2)) ) /
    ( 2.-2.*sqrt(muj2) ) * scale;
  
  return ptmax > 0.*GeV ? ptmax : 0.*GeV;
}

// NOTE: bounds calculated at this step may be too loose
pair<double,double> FFMassiveInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {

  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
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

bool FFMassiveInvertedTildeKinematics::ptzAllowed(pair<Energy,double> ptz) const {
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / lastScale() );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / lastScale() );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / lastScale() );
  
  Energy pt = ptz.first;
  double z = ptz.second;
  
  double y = ( sqr( pt / lastScale() ) + sqr(1.-z)*mui2 + z*z*mu2 ) /
      (z*(1.-z)*(1.-mui2-mu2-muj2));
  
  // check (y,z) phasespace boundary
  // TODO: is y boundary necessary?
  double bar = 1.-mui2-mu2-muj2;
  double ym = 2.*sqrt(mui2)*sqrt(mu2)/bar;
  double yp = 1. - 2.*sqrt(muj2)*(1.-sqrt(muj2))/bar;
  double zm = ( (2.*mui2+bar*y)*(1.-y) - sqrt(y*y-ym*ym)*sqrt(sqr(2.*muj2+bar-bar*y)-4.*muj2) ) /
    ( 2.*(1.-y)*(mui2+mu2+bar*y) );
  double zp = ( (2.*mui2+bar*y)*(1.-y) + sqrt(y*y-ym*ym)*sqrt(sqr(2.*muj2+bar-bar*y)-4.*muj2) ) /
    ( 2.*(1.-y)*(mui2+mu2+bar*y) );
  
  if ( y<ym || y>yp || z<zm || z>zp ) return false;
  return true;
}

pair<Energy,double> FFMassiveInvertedTildeKinematics::generatePtZ(double& jac, const double * r) const {

  double kappaMin = 
    ptCut() != ZERO ?
    sqr(ptCut()/ptMax()) :
    sqr(0.1*GeV/GeV);

  double kappa;

  using namespace RandomHelpers;

  if ( ptCut() > ZERO ) {
    pair<double,double> kw =
      generate(inverse(0.,kappaMin,1.),r[0]);
    kappa = kw.first;
    jac *= kw.second;
  } else {
    pair<double,double> kw =
      generate((piecewise(),
		flat(1e-4,kappaMin),
		match(inverse(0.,kappaMin,1.))),r[0]);
    kappa = kw.first;
    jac *= kw.second;
  }

  Energy pt = sqrt(kappa)*ptMax();

  pair<double,double> zLims = zBounds(pt);

  pair<double,double> zw =
    generate(inverse(0.,zLims.first,zLims.second)+
	     inverse(1.,zLims.first,zLims.second),r[1]);

  double z = zw.first;
  jac *= zw.second;

  jac *= sqr(ptMax()/lastScale());
  
  if( !ptzAllowed(make_pair(pt,z)) ) jac = 0.;

  return make_pair(pt,z);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMassiveInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FFMassiveInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FFMassiveInvertedTildeKinematics::Init() {

  static ClassDocumentation<FFMassiveInvertedTildeKinematics> documentation
    ("FFMassiveInvertedTildeKinematics inverts the final-final tilde "
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMassiveInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigFFMassiveInvertedTildeKinematics("Herwig::FFMassiveInvertedTildeKinematics", "Herwig.so");
