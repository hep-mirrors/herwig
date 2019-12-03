// -*- C++ -*-
//
// FFMassiveInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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

#include "ThePEG/Interface/Switch.h"
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

// Matches Stephen Webster's thesis
bool FFMassiveInvertedTildeKinematics::doMap(const double * r) {

  if ( ptMax() < ptCut() ) {
    jacobian(0.0);
    return false;
  }

  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();

  double mapping = 1.0;
  vector<double> values(6);
  pair<Energy,double> ptz = generatePtZ(mapping, r, &values);
  if ( mapping == 0.0 ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;
  Energy2 pt2 = sqr(pt);
  double z = ptz.second;

  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / lastScale() );
  double muj2 = sqr( realEmissionData()->hardProcessMass() / lastScale() );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / lastScale() );
  double Muij2 = sqr( bornEmitterData()->hardProcessMass() / lastScale() );
  double Muk2 = sqr( bornSpectatorData()->hardProcessMass() / lastScale() );

  // Define the scale
  Energy2 Qijk = sqr(lastScale());

  // Most of the required values have been calculated in ptzAllowed
  double y = values[0];
  double zi = values[1];
  double xk = values[2];
  double xij = values[3];
  double QijN2 = values[4];
  double sijkN = values[5];
  double sijkN2 = sqr(sijkN);

  // Construct reference momenta nk, nij
  Lorentz5Momentum nij = ( sijkN2 / (sijkN2-Muij2*Muk2) )
    * (emitter - (Muij2/sijkN)*spectator);
  Lorentz5Momentum nk = ( sijkN2 / (sijkN2-Muij2*Muk2) )
    * (spectator - (Muk2/sijkN)*emitter);

  // Construct qt
  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum qt = getKt(emitter,spectator,pt,phi);

  // Construct qij, qk, qi and qj
  Lorentz5Momentum qij = xij*nij + (Muij2/(xij*sijkN))*nk;
  Lorentz5Momentum spe = xk*nk + (Muk2/(xk*sijkN))*nij;

  Lorentz5Momentum em = z*qij
    + ((pt2/Qijk + mui2 - z*z*Muij2)/(xij*sijkN*z))*nk + qt;
  Lorentz5Momentum emm = (1.-z)*qij
    + ((pt2/Qijk + muj2 - sqr(1.-z)*Muij2)/(xij*sijkN*(1.-z)))*nk - qt;

  em.setMass(realEmitterData()->hardProcessMass());
  em.rescaleEnergy();
  emm.setMass(realEmissionData()->hardProcessMass());
  emm.rescaleEnergy();
  spe.setMass(realSpectatorData()->hardProcessMass());
  spe.rescaleEnergy();

  // book
  realEmitterMomentum() = em;
  realEmissionMomentum() = emm;
  realSpectatorMomentum() = spe;

  // Calculate the jacobian
  double bar = 1.-mui2-muj2-muk2;
  
  // mapFactor defined as dy dz = mapFactor * dpt2/sqr(lastScale()) dz
  double mapFactor = 0.0;
  mapFactor = y*(sqr(lastScale()) / (pt2 + sqr(1.-z)*mui2*Qijk + sqr(z)*muj2*Qijk))
      * abs(1. - 2.*Muk2*QijN2 / (bar*(1.-y)*xij*xk*sijkN));
  
  // Mapping includes only the variable changes/jacobians
  mapping *= mapFactor;
  jacobian( (Qijk*sqr(bar)/rootOfKallen(1.,Muij2,Muk2)) * mapping
            * (1.-y)/(16.*sqr(Constants::pi)) / sHat() );
  
  // Store the parameters
  subtractionParameters().resize(3);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = zi;
  subtractionParameters()[2] = z;
  
  return true;
}


Energy FFMassiveInvertedTildeKinematics::lastPt() const {
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  
  double y = subtractionParameters()[0];
  double z = subtractionParameters()[2];
  
  Energy ret = scale * sqrt( y * (1.-mui2-muj2-muk2) * z*(1.-z) - sqr(1.-z)*mui2 - sqr(z)*muj2 );

  return ret;
}

double FFMassiveInvertedTildeKinematics::lastZ() const {
  return subtractionParameters()[2];
}    

Energy FFMassiveInvertedTildeKinematics::ptMax() const {
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  
  Energy ptmax = rootOfKallen( mui2, muj2, sqr(1.-sqrt(muk2)) ) /
    ( 2.-2.*sqrt(muk2) ) * scale;
  
  return ptmax > 0.*GeV ? ptmax : 0.*GeV;
}

// NOTE: bounds calculated at this step may be too loose
pair<double,double> FFMassiveInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {

  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
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

// Matches Stephen Webster's thesis
bool FFMassiveInvertedTildeKinematics::ptzAllowed(pair<Energy,double> ptz, vector<double>* values) const {

  Energy pt = ptz.first;
  Energy2 pt2 = sqr(pt);
  double z = ptz.second;

  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / lastScale() );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / lastScale() );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / lastScale() );
  double Muij2 = sqr( bornEmitterData()->hardProcessMass() / lastScale() );
  double Muk2 = sqr( bornSpectatorData()->hardProcessMass() / lastScale() );
  
  // Calculate the scales that we need
  Energy2 Qijk = sqr(lastScale());
  double QijN2 = (pt2/Qijk + (1.-z)*mui2 + z*muj2) / z / (1.-z);
  double sijkN = 0.5*( 1. - Muij2 - Muk2 + sqrt( sqr(1.-Muij2-Muk2) - 4.*Muij2*Muk2 ) );
  double bar = 1.-mui2-muj2-muk2;

  double y = ( pt2/Qijk + sqr(1.-z)*mui2 + z*z*muj2 ) /
      (z*(1.-z)*bar);

  // Calculate the scaling factors, xk and xij
  double lambdaK = 1. + (Muk2/sijkN);
  double lambdaIJ = 1. + (Muij2/sijkN);
  double fac1 = lambdaIJ*lambdaK + (muk2 - QijN2)/sijkN;
  double xk =
    ( fac1 + sqrt( sqr(fac1) - 4.*lambdaIJ*lambdaK*muk2/sijkN ) )
    / 2. / lambdaK ;
  double xij = 1. - muk2*(1.-xk) / xk / sijkN;

  // Calculate zi
  double zi =
    ( z*xij*xk*sijkN + muk2*(pt2/Qijk + mui2) / (z*xij*xk*sijkN) )
    / (1.-y) / bar;
  
  // Limits on zi
  double facA = (2.*mui2+bar*y)/2./(mui2 + muj2 + bar*y);
  double facB =
    sqrt( (sqr(2.*muk2 + bar*(1.-y)) - 4.*muk2) *
          (sqr(bar)*sqr(y) - 4.*mui2*muj2))
    / bar / (1.-y) / (bar*y + 2.*mui2);
  double zim = facA * (1. - facB);
  double zip = facA * (1. + facB);
  
  // check (y,z) phase space boundary
  double ym = 2.*sqrt(mui2)*sqrt(muj2)/bar;
  double yp = 1. - 2.*sqrt(muk2)*(1.-sqrt(muk2))/bar;

  if ( y<ym || y>yp || zi<zim || zi>zip ) return false;

  assert( (*values).size() == 6);
  (*values)[0] = y;
  (*values)[1] = zi;
  (*values)[2] = xk;
  (*values)[3] = xij;
  (*values)[4] = QijN2;
  (*values)[5] = sijkN;

  return true;
}


// This is used to generate pt and z
pair<Energy,double> FFMassiveInvertedTildeKinematics::generatePtZ(double& jac, const double * r, vector<double> * values) const {

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
  
  if( !ptzAllowed(make_pair(pt,z), values )) jac = 0.;

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
     "kinematics involving a massive particle.");
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMassiveInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigFFMassiveInvertedTildeKinematics("Herwig::FFMassiveInvertedTildeKinematics", "Herwig.so");
