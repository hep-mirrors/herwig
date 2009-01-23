// -*- C++ -*-
//
// NLO2to2Kinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLO2to2Kinematics class.
//

#include "NLO2to2Kinematics.h"

using namespace Herwig;

born2to2Kinematics::born2to2Kinematics() {}
    
born2to2Kinematics::born2to2Kinematics(vector<Lorentz5Momentum> Momenta, double xp, double xm) {
    
  // Leading order momentum fractions and associated etabar's:
  xpb_   = xp;
  etapb_ = sqrt(1.-xpb_);
  xmb_   = xm; 
  etamb_ = sqrt(1.-xmb_);
  
  // The Born momenta according to the notation of the FMNR papers,
  // in the diboson centre of mass frame:
  p1b_ = Momenta[0];
  p2b_ = Momenta[1];
  k1b_ = Momenta[2];
  k2b_ = Momenta[3];

  // The diboson invariant mass / shat, that and uhat:
  sb_ = (p1b_+p2b_).m2();
  tb_ = (p1b_-k1b_).m2();
  ub_ = (p1b_-k2b_).m2();
  
  // Masses of the final state bosons:
  k12b_ = k1b_.m2();
  k22b_ = k2b_.m2();

  // Boost the momenta so they are in the diboson centre of mass frame:
  Lorentz5Momentum ktot(k1b_+k2b_);
  Boost CMSBoostb(-ktot.boostVector());
  p1b_.boost(CMSBoostb);
  p2b_.boost(CMSBoostb);
  k1b_.boost(CMSBoostb);
  k2b_.boost(CMSBoostb);
  
  // The diboson rapidity:
  Yb_ = 0.5*log(xpb_/xmb_);
  // Note Yb_ = + lastY() if flipped_ = false 
  // but  Yb_ = - lastY() if flipped_ = true.
  // Yb_ is always defined with the quark travelling in the +z direction!
    
  // Polar and azimuthal angles of the dibosons in their rest frame:
  theta1b_ = acos(k1b_.z()/k1b_.vect().mag());
  theta2b_ = atan2(k1b_.x(),k1b_.y());

}

void born2to2Kinematics::sanityCheck() const { 
  // The masses and angles can be used to reconstruct the momenta in 
  // the diboson centre of mass frame as follows (as a sanity check):
  Energy4 lambda = sb_*sb_+k12b_*k12b_+k22b_*k22b_
                 -2.*sb_*k12b_- 2.*sb_*k22b_ - 2.*k12b_*k22b_;
  Energy p  = sqrt(lambda)/2./sqrt(sb_);
  Energy pt = p*sin(theta1b_); 
  Energy pl = p*cos(theta1b_); 
  Lorentz5Momentum p1b(0.*GeV,0.*GeV, p1b_.t(),p1b_.t(),0.*GeV);
  Lorentz5Momentum p2b(0.*GeV,0.*GeV,-p2b_.t(),p2b_.t(),0.*GeV);
  Lorentz5Momentum k1b( pt*sin(theta2b_), pt*cos(theta2b_), pl,
			sqrt(p*p+k12b_),sqrt(k12b_));
  Lorentz5Momentum k2b(-pt*sin(theta2b_),-pt*cos(theta2b_),-pl,
    		       sqrt(p*p+k22b_),sqrt(k22b_));

  // Checks showed these momenta agreed with p1b_,p2b_,k1b_,k2b_ to within 
  // 0.01 eV. These therefore correspond to theta's of Frixione et al.

  Lorentz5Momentum p1bdiff(p1b_-p1b);
  Lorentz5Momentum p2bdiff(p2b_-p2b);
  Lorentz5Momentum k1bdiff(k1b_-k1b);
  Lorentz5Momentum k2bdiff(k2b_-k2b);

  if(p1bdiff.x()/GeV>1.e-9||p1bdiff.y()/GeV>1.e-9||p1bdiff.z()/GeV>1.e-9||
     p1bdiff.t()/GeV>1.e-9||p1bdiff.tau()/GeV>1.e-9) {
      cout << "\n\n\nborn2to2Kinematics:\n";
      cout << "p1b_    = " << p1b_/GeV    << endl;
      cout << "p1b     = " << p1b/GeV     << endl;
      cout << "p1bdiff = " << p1bdiff/GeV << endl;
  }
  if(p2bdiff.x()/GeV>1.e-9||p2bdiff.y()/GeV>1.e-9||p2bdiff.z()/GeV>1.e-9||
     p2bdiff.t()/GeV>1.e-9||p2bdiff.tau()/GeV>1.e-9) {
      cout << "\n\n\nborn2to2Kinematics:\n";
      cout << "p2b_    = " << p2b_/GeV    << endl;
      cout << "p2b     = " << p2b/GeV     << endl;
      cout << "p2bdiff = " << p2bdiff/GeV << endl;
  }
  if(k1bdiff.x()/GeV>1.e-9||k1bdiff.y()/GeV>1.e-9||k1bdiff.z()/GeV>1.e-9||
     k1bdiff.t()/GeV>1.e-9||k1bdiff.tau()/GeV>1.e-9) {
      cout << "\n\n\nborn2to2Kinematics:\n";
      cout << "k1b_    = " << k1b_/GeV    << endl;
      cout << "k1b     = " << k1b/GeV     << endl;
      cout << "k1bdiff = " << k1bdiff/GeV << endl;
  }
  if(k2bdiff.x()/GeV>1.e-9||k2bdiff.y()/GeV>1.e-9||k2bdiff.z()/GeV>1.e-9||
     k2bdiff.t()/GeV>1.e-9||k2bdiff.tau()/GeV>1.e-9) { 
      cout << "\n\n\nborn2to2Kinematics:\n";
      cout << "k2b_    = " << k2b_/GeV    << endl;
      cout << "k2b     = " << k2b/GeV     << endl;
      cout << "k2bdiff = " << k2bdiff/GeV << endl;
  }

  return;
}


real2to3Kinematics::real2to3Kinematics() {}
    
real2to3Kinematics::real2to3Kinematics(born2to2Kinematics bornVariables, 
				       double xt, double y) {

  // Store the born2to2Kinematics object
  bornVariables_ = bornVariables;

  // Store the `raw' radiative variables
  xt_ = xt;
  y_  = y ;

  // First we get the lower limit on the x integration, xbar:
  if(y_== 1.)      xbar_ = bornVariables.xpb(); 
  else if(y_==-1.) xbar_ = bornVariables.xmb();
  else {
    double xpb(bornVariables.xpb()), xmb(bornVariables.xmb());
    double xp2b(xpb*xpb)           , xm2b(xmb*xmb);
    double omy(1.-y_);
    double opy(1.+y_);
    double xbar1=2.*opy*xp2b/
	(sqrt(sqr(1.+xp2b)*sqr(omy)+16.*y_*xp2b)+omy*(1.-xpb)*(1.+xpb));
    double xbar2=2.*omy*xm2b/
	(sqrt(sqr(1.+xm2b)*sqr(opy)-16.*y_*xm2b)+opy*(1.-xmb)*(1.+xmb));
    xbar_ = max(xbar1,xbar2);
  }

  // Then we calculate x from \tilde{x}:
  xr_ = xt_==1. ? 1. : xbar_+(1.-xbar_)*xt_;

  if(xr_== 1.)      xpr_ = bornVariables.xpb();
  else if(y_ ==-1.) xpr_ = bornVariables.xpb();
  else if(y_ == 1.) xpr_ = bornVariables.xpb()/xr_;
  else xpr_ = (bornVariables.xpb()/sqrt(xr_))
            * sqrt((2.-(1.-xr_)*(1.-y_))/(2.-(1.-xr_)*(1.+y_)));

  if(xr_== 1.)      xmr_ = bornVariables.xmb();
  else if(y_ ==-1.) xmr_ = bornVariables.xmb()/xr_;
  else if(y_ == 1.) xmr_ = bornVariables.xmb();
  else xmr_ = (bornVariables.xmb()/sqrt(xr_))
            * sqrt((2.-(1.-xr_)*(1.+y_))/(2.-(1.-xr_)*(1.-y_)));

  // The diboson invariant mass is preserved as are the individual 
  // diboson masses, as are theta1 and theta2:
  s2r_     = bornVariables.sb();
  k12r_    = bornVariables.k12b();
  k22r_    = bornVariables.k22b();
  theta1r_ = bornVariables.theta1b();
  theta2r_ = bornVariables.theta2b();

  if(xt_==1.) {
    // Now determine the variables of Frixione et al. for this x and y.
    // Eq.2.6 of WZ paper NPB 383(1992) 3-44): 
    sr_  = s2r_;
    tkr_ = 0.*GeV2;
    ukr_ = 0.*GeV2;
    // Eq.2.12 of WZ paper NPB 383(1992) 3-44:     
    cpsir_  =  -1.;
    cpsiprr_= 999.;
  } 
  else if(y_== 1.) {
    // Now determine the variables of Frixione et al. for this x and y.
    // Eq.2.6 of WZ paper NPB 383(1992) 3-44): 
    sr_  = s2r_/xr_;
    tkr_ =  0.*GeV2;
    ukr_ = -sr_*(1.-xr_);
    // Eq.2.12 of WZ paper NPB 383(1992) 3-44: 
    cpsir_  = -1.;
    cpsiprr_=  1.;
  } 
  else if(y_==-1.) { 
    // Now determine the variables of Frixione et al. for this x and y.
    // Eq.2.6 of WZ paper NPB 383(1992) 3-44): 
    sr_  = s2r_/xr_;
    tkr_ = -sr_*(1.-xr_);
    ukr_ =  0.*GeV2;
    // Eq.2.12 of WZ paper NPB 383(1992) 3-44: 
    cpsir_  = -1.;
    cpsiprr_= -1.;
  } 
  else {
    // Now determine the variables of Frixione et al. for this x and y.
    // Eq.2.6 of WZ paper NPB 383(1992) 3-44): 
    sr_  = s2r_/xr_;
    tkr_ = -0.5*sr_*(1.-xr_)*(1.-y_);
    ukr_ = -0.5*sr_*(1.-xr_)*(1.+y_);
    // Eq.2.12 of WZ paper NPB 383(1992) 3-44: 
    cpsir_  =1.-sr_ /2./((sr_+tkr_)/2./sqrt(s2r_))/((  sr_+ukr_)/2./sqrt(s2r_));
    cpsiprr_=1.+tkr_/2./((sr_+tkr_)/2./sqrt(s2r_))/((-tkr_-ukr_)/2./sqrt(s2r_));
  }

  // Remainder of Eq.2.12 of WZ paper NPB 383(1992) 3-44: 
  betaxr_ = sqrt(1.-sqr(sqrt(k12r_)+sqrt(k22r_))/s2r_)
          * sqrt(1.-sqr(sqrt(k12r_)-sqrt(k22r_))/s2r_);
  v1r_    = betaxr_/(1.-(k22r_-k12r_)/s2r_); 
  v2r_    = betaxr_/(1.+(k22r_-k12r_)/s2r_); 
  
  // Eq.2.13 of WZ paper NPB 383(1992) 3-44: 
  if(xt_==1.) {
    q1r_ = k12r_ - 0.5*(sr_     )*betaxr_/v1r_*(1.-v1r_*cos(theta1r_));
    q2r_ = k22r_ - 0.5*(sr_     )*betaxr_/v2r_*(1.-v2r_*cos(theta1r_));
  } 
  else if(y_== 1.) {
    q1r_ = k12r_ - 0.5*(sr_     )*betaxr_/v1r_*(1.-v1r_*cos(theta1r_));
    q2r_ = k22r_ - 0.5*(sr_*xr_ )*betaxr_/v2r_*(1.-v2r_*cos(theta1r_));
  } 
  else if(y_==-1.) { 
    q1r_ = k12r_ - 0.5*(sr_*xr_ )*betaxr_/v1r_*(1.-v1r_*cos(theta1r_));
    q2r_ = k22r_ - 0.5*(sr_     )*betaxr_/v2r_*(1.-v2r_*cos(theta1r_));
  } 
  else {
    q1r_ = k12r_ - 0.5*(sr_+tkr_)*betaxr_/v1r_*(1.-v1r_*cos(theta1r_));
    q2r_ = k22r_ - 0.5*(sr_+ukr_)*betaxr_/v2r_*(1.+v2r_*cos(theta2r_)
		       *sin(theta1r_)*sqrt(1.-cpsir_)*sqrt(1.+cpsir_)
	                                          +v2r_*cos(theta1r_)*cpsir_);
  }

  // Eq.2.7 of WZ paper NPB 383(1992) 3-44 (s2 was already defined above): 
  q1hatr_ = k12r_ + k22r_ - sr_  - tkr_ - q1r_;
  q2hatr_ = k12r_ + k22r_ - sr_  - ukr_ - q2r_;
  w1r_    = k12r_ - q1r_  + q2r_ - tkr_;
  w2r_    = k22r_ + q1r_  - q2r_ - ukr_;
  
  

}

void real2to3Kinematics::sanityCheck() const {

  // Sanity check: reconstruct the momenta and check agreement with 
  // k1b_, k2b_. These should not have changed modulo +/- conventions
  // for things like sin(psi) and sin(psi^prime) (cf Eq.2.12 of WZ paper).
  Energy p10r( (sr_ +tkr_)/2./sqrt(s2r_));
  Energy p20r( (sr_ +ukr_)/2./sqrt(s2r_));
  Energy k0r (-(ukr_+tkr_)/2./sqrt(s2r_));
  double spsir   = sqrt(1-cpsir_  )*sqrt(1.+cpsir_  );
  double spsiprr = sqrt(1-cpsiprr_)*sqrt(1.+cpsiprr_);
  Lorentz5Momentum p1r(0.*GeV,0.*GeV,p10r,p10r,0.*GeV);
  Lorentz5Momentum p2r(0.*GeV, p20r*spsir  ,p20r*cpsir_  ,p20r,0.*GeV);
  Lorentz5Momentum kr (0.*GeV, k0r *spsiprr,k0r *cpsiprr_,k0r ,0.*GeV);
  Lorentz5Momentum k1r(sqrt(s2r_)*v1r_*sin(theta2r_)*sin(theta1r_),
		       sqrt(s2r_)*v1r_*cos(theta2r_)*sin(theta1r_),
		       sqrt(s2r_)*v1r_*cos(theta1r_),
		       sqrt(s2r_),sqrt(k12r_));
  Lorentz5Momentum k2r(sqrt(s2r_)*-v2r_*sin(theta2r_)*sin(theta1r_),
		       sqrt(s2r_)*-v2r_*cos(theta2r_)*sin(theta1r_),
		       sqrt(s2r_)*-v2r_*cos(theta1r_),
		       sqrt(s2r_),sqrt(k22r_));
  k1r *= betaxr_/2./v1r_;
  k2r *= betaxr_/2./v2r_;
  k1r.setTau(sqrt(k12r_));
  k2r.setTau(sqrt(k22r_));

  // Check that everything is on shell (nearly):
  if((p1r.m()-p1r.mass())/MeV>5.||(p2r.m()-p2r.mass())/MeV>5.||
     (kr.m() -kr.mass() )/MeV>5.||
     (k1r.m()-k1r.mass())/MeV>5.||(k2r.m()-k2r.mass())/MeV>5.) {
      cout << "\nreal2to3Kinematics off-shell particle(s) reconstructed:\n";
      cout << "p1r  = " << p1r       /GeV << "   "
	   << "mass = " << p1r.mass()/GeV << "   m = " << p1r.m()/GeV << endl;
      cout << "p2r  = " << p2r       /GeV << "   "
	   << "mass = " << p2r.mass()/GeV << "   m = " << p2r.m()/GeV << endl;
      cout << "kr   = " << kr        /GeV << "   "
	   << "mass = " << kr.mass() /GeV << "   m = " << kr.m() /GeV << endl;
      cout << "k1r  = " << k1r       /GeV << "   "
	   << "mass = " << k1r.mass()/GeV << "   m = " << k1r.m()/GeV << endl;
      cout << "k2r  = " << k2r       /GeV << "   "
	   << "mass = " << k2r.mass()/GeV << "   m = " << k2r.m()/GeV << endl;
  }

  // Check that the total momentum is conserved:
  Lorentz5Momentum total;
  total = p1r+p2r-kr-k1r-k2r;
  if(total.x()/MeV>0.1||total.y()/MeV>0.1||total.z()/MeV>0.1||
     total.t()/MeV>0.1)
     cout << "\nreal2to3Kinematics momentum imbalance = " << total/GeV << endl;

  // OK rescale all spatial components so the invariant length
  // equals exactly the invariant length element of the 5Vector.
  p1r.rescaleRho();
  p2r.rescaleRho();
  kr.rescaleRho();
  k1r.rescaleRho();
  k2r.rescaleRho();

  // Check that the total momentum is still conserved after the rescaling:
  total = p1r+p2r-kr-k1r-k2r;
  if(total.x()/MeV>0.1||total.y()/MeV>0.1||total.z()/MeV>0.1||
     total.t()/MeV>0.1)
     cout << "\nreal2to3Kinematics momentum imbalance = " << total/GeV << endl;

  // Check that the final state momenta k1 and k2 are the same as
  // they are in the 2->2 process:
  Lorentz5Momentum k1diff(bornVariables_.k1b()-k1r);
  Lorentz5Momentum k2diff(bornVariables_.k2b()-k2r);
  if(k1diff.x()/GeV>1.e-9||k1diff.y()/GeV>1.e-9||k1diff.z()/GeV>1.e-9||
     k1diff.t()/GeV>1.e-9||k1diff.tau()/GeV>1.e-9) {
      cout << "\n\n\nreal2to3Kinematics:\n";
      cout << "k1b    = " << bornVariables_.k1b()/GeV << endl;
      cout << "k1r    = " << k1r   /GeV << endl;
      cout << "k1diff = " << k1diff/GeV << endl;
  }
  if(k2diff.x()/GeV>1.e-9||k2diff.y()/GeV>1.e-9||k2diff.z()/GeV>1.e-9||
     k2diff.t()/GeV>1.e-9||k2diff.tau()/GeV>1.e-9) {
      cout << "\n\n\nreal2to3Kinematics:\n";
      cout << "k2b    = " << bornVariables_.k2b()/GeV << endl;
      cout << "k2r    = " << k2r   /GeV << endl;
      cout << "k2diff = " << k2diff/GeV << endl;
  }

  // Check also that for y_=+/-1 you also get Born-like initial state momenta
  if(xt_!=1.&&y_== 1.) {
    Lorentz5Momentum p1bdiff = (p1r-kr-bornVariables_.p1b());
    if(p1bdiff.x()/GeV>1.e-9||p1bdiff.y()/GeV>1.e-9||p1bdiff.z()/GeV>1.e-9||
       p1bdiff.t()/GeV>1.e-9||p1bdiff.tau()/GeV>1.e-9) {
	cout << "\n\n\nreal2to3Kinematics: error for y_= 1.\n";
	cout << "p1r-kr     = " << (p1r-kr)                     /GeV << endl;
	cout << "p1b        = " <<  bornVariables_.p1b()        /GeV << endl;
	cout << "p1r-kr-p1b = " << (p1r-kr-bornVariables_.p1b())/GeV << endl;
    }
    Lorentz5Momentum p2bdiff = (p2r-bornVariables_.p2b());
    if(p2bdiff.x()/GeV>1.e-9||p2bdiff.y()/GeV>1.e-9||p2bdiff.z()/GeV>1.e-9||
       p2bdiff.t()/GeV>1.e-9||p2bdiff.tau()/GeV>1.e-9) {
	cout << "\n\n\nreal2to3Kinematics: error for y_=-1.\n";
	cout << "p2r     = " <<  p2r                      /GeV << endl;
	cout << "p2b     = " <<  bornVariables_.p2b()     /GeV << endl;
	cout << "p2r-p2b = " << (p2r-bornVariables_.p2b())/GeV << endl;
    }
  }
  if(xt_!=1.&&y_==-1.) {
    Lorentz5Momentum p2bdiff = (p2r-kr-bornVariables_.p2b());
    if(p2bdiff.x()/GeV>1.e-9||p2bdiff.y()/GeV>1.e-9||p2bdiff.z()/GeV>1.e-9||
       p2bdiff.t()/GeV>1.e-9||p2bdiff.tau()/GeV>1.e-9) {
	cout << "\n\n\nreal2to3Kinematics: error for y_=-1.\n";
	cout << "p2r-kr     = " << (p2r-kr)                     /GeV << endl;
	cout << "p2b        = " <<  bornVariables_.p2b()        /GeV << endl;
	cout << "p2r-kr-p2b = " << (p2r-kr-bornVariables_.p2b())/GeV << endl;
    }
    Lorentz5Momentum p1bdiff = (p1r-bornVariables_.p1b());
    if(p1bdiff.x()/GeV>1.e-9||p1bdiff.y()/GeV>1.e-9||p1bdiff.z()/GeV>1.e-9||
       p1bdiff.t()/GeV>1.e-9||p1bdiff.tau()/GeV>1.e-9) {
	cout << "\n\n\nreal2to3Kinematics: error for y_=-1.\n";
	cout << "p1r     = " <<  p1r                      /GeV << endl;
	cout << "p1b     = " <<  bornVariables_.p1b()     /GeV << endl;
	cout << "p1r-p1b = " << (p1r-bornVariables_.p1b())/GeV << endl;
    }
  }

  // And that for xt_=1 you also get the Born initial state momenta
  if(xt_== 1.) {
    Lorentz5Momentum p1bdiff = (p1r-bornVariables_.p1b());
    if(p1bdiff.x()/GeV>1.e-9||p1bdiff.y()/GeV>1.e-9||p1bdiff.z()/GeV>1.e-9||
       p1bdiff.t()/GeV>1.e-9||p1bdiff.tau()/GeV>1.e-9) {
	cout << "\n\n\nreal2to3Kinematics: error for xt_= 1.\n";
	cout << "p1r     = " <<  p1r                      /GeV << endl;
	cout << "p1b     = " <<  bornVariables_.p1b()     /GeV << endl;
	cout << "p1r-p1b = " << (p1r-bornVariables_.p1b())/GeV << endl;
    }
    Lorentz5Momentum p2bdiff = (p2r-bornVariables_.p2b());
    if(p2bdiff.x()/GeV>1.e-9||p2bdiff.y()/GeV>1.e-9||p2bdiff.z()/GeV>1.e-9||
       p2bdiff.t()/GeV>1.e-9||p2bdiff.tau()/GeV>1.e-9) {
	cout << "\n\n\nreal2to3Kinematics: error for xt_= 1.\n";
	cout << "p2r     = " <<  p2r                      /GeV << endl;
	cout << "p2b     = " <<  bornVariables_.p2b()     /GeV << endl;
	cout << "p2r-p2b = " << (p2r-bornVariables_.p2b())/GeV << endl;
    }
  }

  return;
}

