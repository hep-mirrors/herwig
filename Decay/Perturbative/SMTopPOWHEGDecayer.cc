// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMTopPOWHEGDecayer class.
//

#include "SMTopPOWHEGDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"

using namespace Herwig;


SMTopPOWHEGDecayer::SMTopPOWHEGDecayer() : mt_(ZERO), w_(0.), b_(0.), w2_(0.),
					   b2_(0.), pTmin_(GeV), pT_(ZERO)
{}

IBPtr SMTopPOWHEGDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SMTopPOWHEGDecayer::fullclone() const {
  return new_ptr(*this);
}


void SMTopPOWHEGDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(pTmin_,GeV);
}

void SMTopPOWHEGDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(pTmin_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMTopPOWHEGDecayer,SMTopDecayer>
describeHerwigSMTopPOWHEGDecayer("Herwig::SMTopPOWHEGDecayer", "HwPerturbativeDecay.so");

void SMTopPOWHEGDecayer::Init() {

  static ClassDocumentation<SMTopPOWHEGDecayer> documentation
    ("There is no documentation for the SMTopPOWHEGDecayer class");

  static Parameter<SMTopPOWHEGDecayer,Energy> interfacepTmin
    ("pTmin",
     "Minimum transverse momentum from gluon radiation",
     &SMTopPOWHEGDecayer::pTmin_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
}


RealEmissionProcessPtr SMTopPOWHEGDecayer::generateHardest(RealEmissionProcessPtr born) {
  PPtr top = born->bornIncoming()[0];
  // get the bottom and W
  assert(born->bornOutgoing().size()==2);
  PPtr bottom = born->bornOutgoing()[0];
  PPtr Wboson = born->bornOutgoing()[1];
  bool order = bottom->id()!=ParticleID::b;
  if(order) swap(bottom,Wboson);
  // masses of the particles
  mt_ = top   ->momentum().mass();
  w_  = Wboson->momentum().mass() / mt_;
  b_  = bottom->momentum().mass() / mt_; 
  w2_ = sqr(w_);
  b2_ = sqr(b_);
  // find rotation fgrom lab to frame with W along -z
  LorentzRotation eventFrame( top->momentum().findBoostToCM() );
  Lorentz5Momentum pspectator = eventFrame*Wboson->momentum();
  eventFrame.rotateZ( -pspectator.phi() );
  eventFrame.rotateY( -pspectator.theta() - Constants::pi );
  //invert it
  eventFrame.invert();
  //generate the hard emission
  vector<Lorentz5Momentum> momenta = hardMomenta();
  // if no emission return
  if(momenta.empty()) {
    born->pT()[ShowerInteraction::QCD] = pTmin_;
    return born;
  }
  // rotate momenta back to the lab
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    momenta[ix] *= eventFrame;  
  }
  // create new Particles
  PPtr newTop = top   ->dataPtr()->produceParticle(top->momentum());
  born->incoming().push_back(newTop);
  PPtr newb = bottom->dataPtr()->produceParticle(momenta[1]);
  PPtr newW = Wboson->dataPtr()->produceParticle(momenta[2]);
  PPtr newg = getParticleData(ParticleID::g)->produceParticle(momenta[3]);
  // colour flow
  newg->incomingColour(newTop,top->id()<0);
  newg->colourConnect (newb  ,top->id()<0);
  // outgoing particles
  if(!order) {
    born->outgoing().push_back(newb);
    born->outgoing().push_back(newW);
    born->emitter(1);
  }
  else {
    born->outgoing().push_back(newW);
    born->outgoing().push_back(newb);
    born->emitter(2);
  }
  born->outgoing().push_back(newg);
  // boost for the W
  LorentzRotation trans(Wboson->momentum().findBoostToCM());
  trans.boost(momenta[2].boostVector());
  born->transformation(trans);
  born->spectator(0);
  born->emitted(3);
  born->pT()[ShowerInteraction::QCD] = pT_;
  born->interaction(ShowerInteraction::QCD);
  return born;
}

vector<Lorentz5Momentum>  SMTopPOWHEGDecayer::hardMomenta() {
  double C    = 6.3;
  double ymax = 10.;
  double ymin = -ymax;
  
  vector<Lorentz5Momentum> particleMomenta (4);
  Energy2 lambda = sqr(mt_)* sqrt( 1. + sqr(w2_) + sqr(b2_) - 2.*w2_ - 2.*b2_ - 2.*w2_*b2_);  

  //Calculate A
  double A = (ymax - ymin) * C * (coupling()->overestimateValue() / (2.*Constants::pi));
   
  Energy pTmax = mt_* (sqr(1.-w_) - b2_) / (2.*(1.-w_));
  if (pTmax < pTmin_) particleMomenta.clear();

  while (pTmax >= pTmin_) {  
    //Generate pT, y and phi values
    Energy pT = pTmax * pow(UseRandom::rnd() , (1./A));  
    
    if (pT < pTmin_) {particleMomenta.clear(); break;}

    double phi = UseRandom::rnd() * Constants::twopi;
    double y   = ymin + UseRandom::rnd() * (ymax-ymin);

    double weight[2] = {0.,0.};
    double xw[2], xb[2], xb_z[2], xg;
    
    for (unsigned int j=0; j<2; j++) {
      //Check if the momenta are physical
      bool physical = calcMomenta(j, pT, y, phi, xg, xw[j], xb[j], xb_z[j], 
				  particleMomenta);
      if (! physical) continue;
      
      //Check if point lies within phase space
      bool inPS = psCheck(xg, xw[j]);
      if (! inPS) continue;
      
      //Calculate the ratio R/B
      double meRatio = matrixElementRatio(particleMomenta);
      
      //Calculate jacobian
      Energy2 denom = (mt_ - particleMomenta[3].e()) * 
	               particleMomenta[2].vect().mag() -
		       particleMomenta[2].e() * particleMomenta[3].z(); 

      InvEnergy2 J  = (particleMomenta[2].vect().mag2()) / (2.* lambda * denom);
      //Calculate weight
      weight[j] = meRatio * fabs(sqr(pT)*J) * coupling()->ratio(pT*pT) / C;    
     }

    //Accept point if weight > R
    if (weight[0] + weight[1] > UseRandom::rnd()) {
      if (weight[0] > (weight[0] + weight[1])*UseRandom::rnd()) {
	particleMomenta[1].setE( (mt_/2.)*xb  [0]);
	particleMomenta[1].setZ( (mt_/2.)*xb_z[0]);
	particleMomenta[2].setE( (mt_/2.)*xw  [0]);
	particleMomenta[2].setZ(-(mt_/2.)*sqrt(sqr(xw[0])-4.*w2_));
      }
      else {
	particleMomenta[1].setE( (mt_/2.)*xb  [1]);
	particleMomenta[1].setZ( (mt_/2.)*xb_z[1]);
	particleMomenta[2].setE( (mt_/2.)*xw  [1]);
	particleMomenta[2].setZ(-(mt_/2.)*sqrt(sqr(xw[1])-4.*w2_));
      }
      pT_ = pT;
      break;   
    }
    //If there's no splitting lower the pT
    pTmax = pT; 
  }
  return particleMomenta;
}

bool SMTopPOWHEGDecayer::deadZoneCheck(double xw, double xg){

  //veto events not in the dead cone 
  double Lambda = sqrt(1. + sqr(w2_) + sqr(b2_) - 2.*w2_ - 2.*b2_ - 2.*w2_*b2_);
  double kappa = b2_ + 0.5*(1. - w2_ + b2_ + Lambda);
  //invert xw for z values
  double A =  1.;
  double B = -1.;
  double C =  (1.+w2_-b2_-xw)/kappa;
  if((sqr(B) - 4.*A*C) >= 0.){
    double z[2];
    z[0] = (-B + sqrt(sqr(B) - 4.*A*C))/(2.*A);
    z[1] = (-B - sqrt(sqr(B) - 4.*A*C))/(2.*A);
    double r = 0.5*(1. + b2_/(1. + w2_- xw));
    double xg_lims [2];
    xg_lims[0] = (2. - xw)*(1.-r) - (z[0]-r)*sqrt(sqr(xw) - 4.*w2_);
    xg_lims[1] = (2. - xw)*(1.-r) - (z[1]-r)*sqrt(sqr(xw) - 4.*w2_);
    double xg_low_lim = min(xg_lims[0], xg_lims[1]);
    double xg_upp_lim = max(xg_lims[0], xg_lims[1]);
    if (xg>=xg_low_lim && xg<=xg_upp_lim) return false;
  }

  double kappa_t = 1. + 0.5*(1. - w2_ + b2_ + Lambda);
  double z = 1. - xg/kappa_t; 
  double u = 1. + w2_ - b2_ - (1.-z)*kappa_t;
  double y = 1. - (1.-z)*(kappa_t-1.);
  if (sqr(u) - 4.*w2_*y*z >= 0.){
    double v = sqrt(sqr(u) - 4.*w2_*y*z);
    double xw_lim = (u + v) / (2.*y) + (u - v) / (2.*z);
    if (xw <= xw_lim) return false;
  }
  else if (sqr(u) - 4.*w2_*y*z < 0.){
    double xg_lim = (8.*w2_ -2.*xw*(1-b2_+w2_))/(4.*w2_-2.*xw);
    if (xg>=xg_lim) return false;
  }

  return true;
}


double SMTopPOWHEGDecayer::matrixElementRatio( 
				    vector<Lorentz5Momentum> particleMomenta) {
  
  double f  = (1. + sqr(b2_) - 2.*sqr(w2_) + w2_ + w2_*b2_ - 2.*b2_);
  double Nc = standardModel()->Nc();
  double Cf = (sqr(Nc) - 1.) / (2.*Nc);  
  double B  = f/w2_;  

  Energy2 PbPg = particleMomenta[1]*particleMomenta[3];
  Energy2 PtPg = particleMomenta[0]*particleMomenta[3];
  Energy2 PtPb = particleMomenta[0]*particleMomenta[1];

  double R = Cf *((-4.*sqr(mt_)*f/w2_) * ((sqr(mt_)*b2_/sqr(PbPg)) + 
		  (sqr(mt_)/sqr(PtPg)) - 2.*(PtPb/(PtPg*PbPg))) +
		  (16. + 8./w2_ + 8.*b2_/w2_) * ((PtPg/PbPg) + (PbPg/PtPg)) -
		  (16./w2_) * (1. + b2_)); 

  return R/B;
}


bool SMTopPOWHEGDecayer::calcMomenta(int j, Energy pT, double y, double phi,
   		        double& xg, double& xw, double& xb, double& xb_z,
			vector<Lorentz5Momentum>& particleMomenta){
  
  //Calculate xg
  xg = 2.*pT*cosh(y) / mt_;
  if (xg>(1. - sqr(b_ + w_)) || xg<0.) return false;

  //Calculate xw
  double xT  = 2.*pT / mt_;
  double A   = 4. - 4.*xg + sqr(xT);
  double B   = 4.*(3.*xg - 2. + 2.*b2_ - 2.*w2_ - sqr(xg) - xg*b2_ + xg*w2_);
  double L   = 1. + sqr(w2_) + sqr(b2_) - 2.*w2_ - 2.*b2_ - 2.*w2_*b2_;
  double det = 16.*( -L*sqr(xT) + pow(xT,4)*w2_ + 2.*xg*sqr(xT)*(1.-w2_-b2_) + 
		      L*sqr(xg) - sqr(xg*xT)*(1. + w2_) + pow(xg,4) + 
		      2.*pow(xg,3)*(- 1. + w2_ + b2_) );

  if (det<0.) return false;
  if (j==0) xw = (-B + sqrt(det))/(2.*A);
  if (j==1) xw = (-B - sqrt(det))/(2.*A);  
  if (xw>(1. + w2_ - b2_) || xw<2.*w_) return false;
  
  //Calculate xb
  xb = 2. - xw - xg;     
  if (xb>(1. + b2_ - w2_) || xb<2.*b_) return false;       

  //Calculate xb_z  
  // Due to precision limitations, it is possible for
  // sqr(xb)-4.*b2_-s to be very small
  // and negative, set to 0 in this case
  double xb_z_mod2 = std::max(sqr(xb) - 4.*b2_  - sqr(xT), 0.);

  double epsilon_p =  -sqrt(sqr(xw) - 4.*w2_) + xT*sinh(y) +
                       sqrt(xb_z_mod2);
  double epsilon_m =  -sqrt(sqr(xw) - 4.*w2_) + xT*sinh(y) - 
                       sqrt(xb_z_mod2);

  if (fabs(epsilon_p) < 1.e-10){
    xb_z =  sqrt(xb_z_mod2);
  }
  else if (fabs(epsilon_m) < 1.e-10){
    xb_z = -sqrt(xb_z_mod2);
  }
  else return false;

  //Check b is on shell
  if (fabs((sqr(xb) - sqr(xT) - sqr(xb_z) - 4.*b2_))>1.e-10) return false;

  //Calculate 4 momenta
  particleMomenta[0].setE   ( mt_);
  particleMomenta[0].setX   ( ZERO);
  particleMomenta[0].setY   ( ZERO);
  particleMomenta[0].setZ   ( ZERO);
  particleMomenta[0].setMass( mt_);

  particleMomenta[1].setE   ( mt_*xb/2.);
  particleMomenta[1].setX   (-pT*cos(phi));
  particleMomenta[1].setY   (-pT*sin(phi));
  particleMomenta[1].setZ   ( mt_*xb_z/2.);
  particleMomenta[1].setMass( mt_*b_);

  particleMomenta[2].setE   ( mt_*xw/2.);
  particleMomenta[2].setX   ( ZERO);
  particleMomenta[2].setY   ( ZERO);
  particleMomenta[2].setZ   (-mt_*sqrt(sqr(xw) - 4.*w2_)/2.);
  particleMomenta[2].setMass( mt_*w_);

  particleMomenta[3].setE   ( pT*cosh(y));
  particleMomenta[3].setX   ( pT*cos(phi));
  particleMomenta[3].setY   ( pT*sin(phi));
  particleMomenta[3].setZ   ( pT*sinh(y));
  particleMomenta[3].setMass( ZERO);
 
  return true;
}


bool SMTopPOWHEGDecayer::psCheck(double xg, double xw) {
  
  //Check is point is in allowed region of phase space
  double xb_star = (1. - w2_ + b2_ - xg) / sqrt(1. - xg);
  double xg_star = xg / sqrt(1. - xg);

  if ((sqr(xb_star) - 4.*b2_) < 1e-10) return false;
  double xw_max = (4. + 4.*w2_ - sqr(xb_star + xg_star) + 
		   sqr(sqrt(sqr(xb_star) - 4.*b2_) + xg_star)) / 4.;
  double xw_min = (4. + 4.*w2_ - sqr(xb_star + xg_star) + 
		   sqr(sqrt(sqr(xb_star) - 4.*b2_) - xg_star)) / 4.;

  if (xw < xw_min || xw > xw_max) return false;

  return true;
}

