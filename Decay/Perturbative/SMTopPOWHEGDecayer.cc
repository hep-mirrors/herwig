// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMTopPOWHEGDecayer class.
//

#include "SMTopPOWHEGDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;


SMTopPOWHEGDecayer::SMTopPOWHEGDecayer() {}

IBPtr SMTopPOWHEGDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SMTopPOWHEGDecayer::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void SMTopPOWHEGDecayer::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void SMTopPOWHEGDecayer::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<SMTopPOWHEGDecayer,SMTopDecayer>
  describeHerwigSMTopPOWHEGDecayer("Herwig::SMTopPOWHEGDecayer", "HwPertrubativeDecay.so");

void SMTopPOWHEGDecayer::Init() {

  static ClassDocumentation<SMTopPOWHEGDecayer> documentation
    ("There is no documentation for the SMTopPOWHEGDecayer class");

}


HardTreePtr SMTopPOWHEGDecayer::generateHardest(ShowerTreePtr) {

  remove("weights.top");
  unsigned int npoint=1000000;
  ofstream file("dalitz.top");
 
  for(unsigned int ix=0; ix < npoint; ++ix) {
    vector<Lorentz5Momentum> momenta = hardMomenta();
    if (momenta.size()==4){
      double xg = 2.*momenta[3].e()/momenta[0].mass();
      double xw = 2.*momenta[2].e()/momenta[0].mass();
      file << xg << "\t" << xw << "\n";
    }
  }
  file.close();
  exit(1);
}


vector<Lorentz5Momentum>  SMTopPOWHEGDecayer::hardMomenta() {

  Energy mt = getParticleData(ParticleID::t)->mass();
  double w  = getParticleData(ParticleID::Wplus)->mass() / mt;
  double b  = getParticleData(ParticleID::b)->mass() / mt;
  
  Energy pTmin = 1.*GeV;

  vector<Lorentz5Momentum> particleMomenta (4);
  double ymin = -10.;
  double ymax = 10.;
  double C = 1.;
  Energy2 lambda = sqr(mt)* sqrt( 1. + pow(w,4) + pow(b,4) - 
			         2.*sqr(w) - 2.*sqr(b) - 2.*sqr(w*b));    

  //Calculate A
  double A = (ymax - ymin) * C * (coupling()->overestimateValue() / 
				 (2.*Constants::pi));
  Energy pTmax = mt* (sqr(1.-w) - sqr(b)) / (2.*(1.-w));

  while (pTmax >= pTmin){  
    //Generate pT, y and phi values
    Energy pT = pTmax * pow(UseRandom::rnd() , (1./A));  
    if (pT < pTmin) {
      particleMomenta.clear(); 
      break;
    }
    double phi = UseRandom::rnd() * Constants::twopi;
    double y = ymin + UseRandom::rnd() * (ymax-ymin);
    
    double weight[2] = {0.,0.};
    double xw[2], xg;
    
    for (unsigned int j=0; j<2; j++) {
      //Check if the momenta are physical
      bool physical = calcMomenta(j, pT, y, phi, xg, xw[j], particleMomenta);
      if (not physical) continue;
      
      //Check if point lies within phase space
      bool inPS = psCheck(xg, xw[j]);
      if (not inPS) continue;
      
      //Calculate the ratio R/B
      Energy2 meRatio = matrixElementRatio (particleMomenta);
      
      //Calculate jacobian 
      double J = (sqr(mt) * particleMomenta[2].vect().mag2()) / 
	         (8. * pow(Constants::pi,3) * lambda * 
		 (particleMomenta[2].vect().mag()*(mt-particleMomenta[3].e()) -
		  particleMomenta[2].e()*particleMomenta[3].z()));
      
      //Calculate weight
      weight[j] = (sqr(pT)/meRatio) * fabs(J) * coupling()->ratio(pT*pT) / C;
    }

    ofstream weights;
    if (weight[0] + weight[1] > 1.){
      weights.open("weights.top", ios::app);
      weights << weight[0]+weight[1] << endl;
    }

    //Accept point if weight > R
    if (weight[0] + weight[1] > UseRandom::rnd()) {
      if (weight[0] > (weight[0] + weight[1])*UseRandom::rnd()) {
	particleMomenta[2].setE((mt/2.)*xw[0]);
	particleMomenta[2].setZ(-(mt/2.)*sqrt(sqr(xw[0])-4.*sqr(w)));
      }
      else {
	particleMomenta[2].setE((mt/2.)*xw[1]);
	particleMomenta[2].setZ(-(mt/2.)*sqrt(sqr(xw[1])-4.*sqr(w)));
      }
      break;   
    }
    //If there's no splitting lower the pT
    pTmax = pT; 
  }
  return particleMomenta;
}


Energy2 SMTopPOWHEGDecayer::matrixElementRatio(
			   vector<Lorentz5Momentum> particleMomenta){

  Energy mt = getParticleData(ParticleID::t)->mass();
  double w  = getParticleData(ParticleID::Wplus)->mass() / mt;
  double b  = getParticleData(ParticleID::b)->mass() / mt;
             
  Energy2 f = sqr(mt) * (1. + pow(b,4) - 2.*pow(w,4) + 
			sqr(w) + sqr(w*b) - 2.*sqr(b));
  double Nc = standardModel()->Nc();
  double Cf = (sqr(Nc) - 1.) / (2.*Nc);  
  Energy2 B = (1./(2.*sqr(w)))*f;  
  double Norm = sqr(Constants::pi)*Cf/2.;

  Energy2 PbPg = particleMomenta[1]*particleMomenta[3];
  Energy2 PtPg = particleMomenta[0]*particleMomenta[3];
  Energy2 PtPb = particleMomenta[0]*particleMomenta[1];

  double R = Norm * ( (-4.*f/sqr(w)) * 
          ((sqr(mt*b)/sqr(PbPg)) + (sqr(mt)/sqr(PtPg)) -2.*(PtPb/(PbPg*PtPg))) +
          (16. + 8.*sqr(1./w) + 8.*sqr(b/w))* ((PtPg/PbPg) + (PbPg/PtPg)) -
	  (16./sqr(w)) * (1. + sqr(b)) ); 
   
  return B/R;
}


bool SMTopPOWHEGDecayer::calcMomenta(int j, Energy pT, double y, double phi,
				     double& xg, double& xw, 
				     vector<Lorentz5Momentum>& particleMomenta){

  Energy mt = getParticleData(ParticleID::t)->mass();
  double w  = getParticleData(ParticleID::Wplus)->mass() / mt;
  double b  = getParticleData(ParticleID::b)->mass() / mt;

  //Calculate xg
  xg = 2.*pT*cosh(y) / mt;
  if (xg>(1. - sqr(b + w)) || xg<0.) return false;

  //Calculate xw
  double xT = 2.*pT / mt;
  double A = 4. - 4.*xg + sqr(xT);
  double B = 4.*(3.*xg - 2. + 2.*sqr(b) - 2.*sqr(w) - sqr(xg) - 
		 xg*sqr(b) + xg*sqr(w));
  double L = 1. + pow(w,4) + pow(b,4) - 2.*sqr(w) - 2.*sqr(b) - 2.*sqr(w*b);
  double det = 16.*(-L*sqr(xT) + pow(xT,4)*sqr(w) + 
		    2.*xg*sqr(xT)*(1. - sqr(w) - sqr(b)) + 
		    L*sqr(xg) - sqr(xg*xT)*(1. + sqr(w)) + 
		    2.*pow(xg,3)*(- 1. + sqr(w) + sqr(b)) + pow(xg,4));

  if (det<0.) return false;
  if (j==0) xw = (-B + sqrt(det))/(2.*A);
  if (j==1) xw = (-B - sqrt(det))/(2.*A);  
  if (xw>(1. + sqr(w) - sqr(b)) || xw<0.) return false;

  //Calculate xb
  double xb = 2. - xw - xg;     
  if (xb>(1. + sqr(b) - sqr(w)) || xb<0.) return false;       

  //Calculate xb_z  
  double xb_z;
  double epsilon_p =  -sqrt(sqr(xw) - 4.*sqr(w)) + xT*sinh(y) +
                       sqrt(sqr(xb) - 4.*sqr(b) - sqr(xT));
  double epsilon_m =  -sqrt(sqr(xw) - 4.*sqr(w)) + xT*sinh(y) - 
                       sqrt(sqr(xb) - 4.*sqr(b) - sqr(xT));

  if (fabs(epsilon_p) < 1.e-6){
    xb_z = sqrt(sqr(xb) - 4.*sqr(b) - sqr(xT));
  }
  else if (fabs(epsilon_m) < 1.e-6){
    xb_z = -sqrt(sqr(xb) - 4.*sqr(b) - sqr(xT));
  }
  else return false;

  //Check b is on shell
  if (fabs((sqr(xb) - sqr(xT) - sqr(xb_z) - 4.*sqr(b)))>1.e-6) return false;

  //Calculate 4 momenta
  particleMomenta[0].setE(mt);
  particleMomenta[0].setMass(mt);

  particleMomenta[1].setE((mt/2.)*xb);
  particleMomenta[1].setX(-pT*cos(phi));
  particleMomenta[1].setY(-pT*sin(phi));
  particleMomenta[1].setZ((mt/2.)*xb_z);
  particleMomenta[1].setMass(mt*b);

  particleMomenta[2].setE((mt/2.)*xw);
  particleMomenta[2].setZ(-(mt/2.)*sqrt(sqr(xw) - 4.*sqr(w)));
  particleMomenta[2].setMass(mt*w);

  particleMomenta[3].setE(pT*cosh(y));
  particleMomenta[3].setX(pT*cos(phi));
  particleMomenta[3].setY(pT*sin(phi));
  particleMomenta[3].setZ(pT*sinh(y));
 
  return true;
}


bool SMTopPOWHEGDecayer::psCheck(double xg, double xw){
  
  Energy mt = getParticleData(ParticleID::t)->mass();
  double w  = getParticleData(ParticleID::Wplus)->mass() / mt;
  double b  = getParticleData(ParticleID::b)->mass() / mt;
            
  //Check is point is in allowed region of phase space
  double xb_star = (1. - sqr(w) + sqr(b) - xg) / sqrt(1. - xg);
  double xg_star = xg / sqrt(1. - xg);

  if ((sqr(xb_star) - 4.*sqr(b)) < 0) return false;
  double xw_max = (4. + 4.*sqr(w) - sqr(xb_star + xg_star) + 
		   sqr(sqrt(sqr(xb_star) - 4.*sqr(b)) + xg_star)) / 4.;
  double xw_min = (4. + 4.*sqr(w) - sqr(xb_star + xg_star) + 
		   sqr(sqrt(sqr(xb_star) - 4.*sqr(b)) - xg_star)) / 4.;

  if (xw < xw_min || xw > xw_max) return false;
  return true;
}

