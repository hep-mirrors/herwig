// -*- C++ -*-
//
// TwoToThreePhaseSpace.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TwoToThreePhaseSpace_H
#define HERWIG_TwoToThreePhaseSpace_H
//
// This is the declaration of the Isospin namespace and values.
//
namespace Herwig {
namespace TwoToThreePhaseSpace {

  /**
   * Generate the \f$2\to3\f$ phase space using a simple mapping for the final-state
   * The returned Jacobian contains all the factors required to give the cross section
   * when multiplied by the spin/colour averaged matrix element
   */
  double twoToThreeFS(Energy ecm, vector<Energy> masses, const double * r,
		      Lorentz5Momentum & q1, Lorentz5Momentum & q2, Lorentz5Momentum & q3,
		      double power=1.) {
    power=0.;
    assert(masses.size()==3);
    double jacobian = 1.;
    unsigned int ioff   = UseRandom::rnd()<0.5 ? 0 : 1;
    unsigned int ispect = ioff==0 ? 1 : 0; 
    // limits for virtual particle mass
    Energy2 mmin(sqr(masses[ioff]+masses[2])),mmax(sqr(ecm-masses[ispect]));
    double rhomin,rhomax;
    if(power==0.) {
      rhomin = mmin/sqr(masses[ioff]);
      rhomax = mmax/sqr(masses[ioff]);
    }
    else if(power==1.) {
      rhomax = log((mmax-sqr(masses[ioff]))/sqr(masses[ioff]));
      rhomin = log((mmin-sqr(masses[ioff]))/sqr(masses[ioff]));
    }
    else {
      rhomin = pow((mmax-sqr(masses[ioff]))/sqr(masses[ioff]),1.-power);
      rhomax = pow((mmin-sqr(masses[ioff]))/sqr(masses[ioff]),1.-power);
      jacobian /= (power-1.);
    }
    double rho = rhomin+r[1]*(rhomax-rhomin);
    Energy2 moff2;
    if(power==0) 
      moff2 = sqr(masses[ioff])*rho;
    else if(power==1)
      moff2 = sqr(masses[ioff])*(exp(rho)+1.);
    else
      moff2 = sqr(masses[ioff])*(pow(rho,1./(1.-power))+1.);
    Energy moff = sqrt(moff2);
    Energy p1,p2;
    try {
      p1 = SimplePhaseSpace::getMagnitude(sqr(ecm), moff, masses[ispect]);
      p2 = SimplePhaseSpace::getMagnitude(moff2,masses[ioff],masses[2]);
    }
    catch ( ImpossibleKinematics & e ) {
      return -1.;
    }
    double cos1 = -1.+2.*r[0];
    double sin1(sqrt(1.-sqr(cos1)));
    double phi1 = Constants::twopi*UseRandom::rnd();
    Lorentz5Momentum poff(sin1*p1*cos(phi1),sin1*p1*sin(phi1),cos1*p1,sqrt(sqr(p1)+moff2),moff);
    q2.setVect(Momentum3(-sin1*p1*cos(phi1),-sin1*p1*sin(phi1),-cos1*p1));
    q2.setMass(masses[ispect]);
    q2.rescaleEnergy();
    bool test=Kinematics::twoBodyDecay(poff,masses[ioff],masses[2],-1.+2*r[2],r[3]*Constants::twopi,q1,q3);
    if(!test) return -1.;
    Energy2 mother2=(q2+q3).m2();
    double D = 2./(pow(sqr(masses[ioff])/(moff2-sqr(masses[ioff])       ),power)+
		   pow(sqr(masses[ispect])/(mother2-sqr(masses[ispect])),power));
    // calculate jacobian
    jacobian *= ecm/moff*(rhomax-rhomin)*sqr(masses[ioff]/ecm)*D*0.125*p1*p2/sqr(ecm)/pow(Constants::twopi,3);
    if(ioff==1) swap(q1,q2);
    return jacobian;
  }
}
}

#endif /* HERWIG_TwoToThreePhaseSpace_H */
