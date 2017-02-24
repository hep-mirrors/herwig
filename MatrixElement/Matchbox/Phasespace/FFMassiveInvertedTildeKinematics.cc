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

#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

using namespace Herwig;

FFMassiveInvertedTildeKinematics::FFMassiveInvertedTildeKinematics() 
  : theFullMapping(false) {}

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
  vector<double> values;
    pair<Energy,double> ptz = generatePtZ(mapping,r, &values);
  if ( mapping == 0.0 ) {
    jacobian(0.0);
    return false;
  }

  // pt and zPrime = qi.nk / (qi+qj).nk are the generated variables
  Energy pt = ptz.first;
  Energy2 pt2 = sqr(pt);
  double zPrime = ptz.second;
  

  // Define the scale
  Energy scale = (emitter+spectator).m();
  Energy2 Qijk = sqr(scale);

  // Most of the required values have been calculated in ptzAllowed
  double y = values[0];
  double z = values[1];
  double xk = values[3];
  double xij = values[4];
  double suijk = values[5];
  double suijk2 = sqr(values[5]);

  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  double Mui2 = sqr( bornEmitterData()->hardProcessMass() / scale );
  double Muj2 = sqr( bornSpectatorData()->hardProcessMass() / scale );


  // Construct reference momenta nk, nij, nt
  Lorentz5Momentum nij = ( suijk2 / (suijk2-Mui2*Muj2) ) * (emitter - (Mui2/suijk)*spectator);
  Lorentz5Momentum nk = ( suijk2 / (suijk2-Mui2*Muj2) ) * (spectator - (Muj2/suijk)*emitter);

  // Following notation in notes, qt = sqrt(wt)*nt
  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum qt = getKt(nij,nk,pt,phi);

  // Construct qij, qk, qi and qj
  Lorentz5Momentum qij = xij*nij + (Mui2/(xij*suijk))*nk;
  Lorentz5Momentum spe = xk*nk + (Muj2/(xk*suijk))*nij;

  Lorentz5Momentum em = zPrime*qij + ((pt2/Qijk + mui2 - zPrime*zPrime*Mui2)/(xij*suijk*zPrime))*nk + qt;
  Lorentz5Momentum emm = (1.-zPrime)*qij + ((pt2/Qijk + mu2 - sqr(1.-zPrime)*Mui2)/(xij*suijk*(1.-zPrime)))*nk - qt;

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
  
  // Store the jacobian
  // The jacobian here corresponds to dpt2 / sqr(lastscale) NOT dpt2 / pt2.
  // This jacobian is the one-particle phase space  
  if ( true ) {
    double bar = 1.-mui2-mu2-muj2;
    mapping *= (1.-y)/(z*(1.-z));
    jacobian( mapping*(sqr(lastScale())/sHat())*bar/rootOfKallen(1.,Mui2,Muj2) / (16.*sqr(Constants::pi)) );
  }

  // Todo: Switch for the full jacobian including the z->zPrime jacobian
    else if ( theFullMapping ) {
      assert ( false );
    }
    
// Store the parameters
  subtractionParameters().resize(3);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = z;
  subtractionParameters()[2] = zPrime;
   
  return true;

}


Energy FFMassiveInvertedTildeKinematics::lastPt() const {
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  
  double y = subtractionParameters()[0];
  double zPrime = subtractionParameters()[2];
  
  Energy ret = scale * sqrt( y * (1.-mui2-mu2-muj2) * zPrime*(1.-zPrime) - sqr(1.-zPrime)*mui2 - sqr(zPrime)*mu2 );

  return ret;
}

double FFMassiveInvertedTildeKinematics::lastZ() const {
  return subtractionParameters()[2];
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
// These apply to zPrime, which is stored in lastZ()
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

bool FFMassiveInvertedTildeKinematics::ptzAllowed(pair<Energy,double> ptz, vector<double>* values) const {

  Energy pt = ptz.first;
  Energy2 pt2 = sqr(pt);
  double zPrime = ptz.second;

  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / lastScale() );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / lastScale() );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / lastScale() );
  double Mui2 = sqr( bornEmitterData()->hardProcessMass() / lastScale() );
  double Muj2 = sqr( bornSpectatorData()->hardProcessMass() / lastScale() );
  
  // Calculate the scales that we need
  Energy2 Qijk = sqr(lastScale());
  double suijk = 0.5*( 1. - Mui2 - Muj2 + sqrt( sqr(1.-Mui2-Muj2) - 4.*Mui2*Muj2 ) );
  double bar = 1.-mui2-mu2-muj2;

  double y = ( pt2/Qijk + sqr(1.-zPrime)*mui2 + zPrime*zPrime*mu2 ) /
      (zPrime*(1.-zPrime)*bar);

  // Calculate A:=xij*w
  double A = (1./(suijk*zPrime*(1.-zPrime))) * ( pt2/Qijk + zPrime*mu2 + (1.-zPrime)*mui2 - zPrime*(1.-zPrime)*Mui2 );

  // Calculate the scaling factors, xk and xij
  double lambdaK = 1. + (Muj2/suijk);
  double lambdaIJ = 1. + (Mui2/suijk);
  double xk = (1./(2.*lambdaK)) * ( (lambdaK + (Muj2/suijk)*lambdaIJ - A) + sqrt( sqr(lambdaK + (Muj2/suijk)*lambdaIJ - A) - 4.*lambdaK*lambdaIJ*Muj2/suijk) );
  double xij = 1. - ( (Muj2/suijk) * (1.-xk) / xk );

  // Calculate z
  double z = ( (zPrime*xij*xk*suijk/2.) + (Muj2/ ( 2.*xk*xij*suijk*zPrime))*(pt2/Qijk + mui2) ) /
    ( (xij*xk*suijk/2.) + (Muj2/(2.*xk*xij))*(Mui2/suijk + A) );
  
  // check (y,z) phasespace boundary
  // TODO: is y boundary necessary?
  double ym = 2.*sqrt(mui2)*sqrt(mu2)/bar;
  double yp = 1. - 2.*sqrt(muj2)*(1.-sqrt(muj2))/bar;
  // These limits apply to z, not zPrime
  double zm = ( (2.*mui2+bar*y)*(1.-y) - sqrt(y*y-ym*ym)*sqrt(sqr(2.*muj2+bar-bar*y)-4.*muj2) ) /
    ( 2.*(1.-y)*(mui2+mu2+bar*y) );
  double zp = ( (2.*mui2+bar*y)*(1.-y) + sqrt(y*y-ym*ym)*sqrt(sqr(2.*muj2+bar-bar*y)-4.*muj2) ) /
    ( 2.*(1.-y)*(mui2+mu2+bar*y) );
  
  if ( y<ym || y>yp || z<zm || z>zp ) return false;

  values->push_back(y);
  values->push_back(z);
  values->push_back(A);
  values->push_back(xk);
  values->push_back(xij);
  values->push_back(suijk);
      
  return true;
}


// This is used to generate pt and zPrime
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
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMassiveInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigFFMassiveInvertedTildeKinematics("Herwig::FFMassiveInvertedTildeKinematics", "Herwig.so");
