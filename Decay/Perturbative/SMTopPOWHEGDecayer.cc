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

using namespace Herwig;


SMTopPOWHEGDecayer::SMTopPOWHEGDecayer() : mt_(ZERO), w_(0.), b_(0.), pTmin_(GeV), pT_(ZERO)
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


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
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


HardTreePtr SMTopPOWHEGDecayer::generateHardest(ShowerTreePtr tree) {
  // get the bottom and W
  assert(tree->outgoingLines().size()==2);
  ShowerProgenitorPtr 
    bProgenitor = tree->outgoingLines(). begin()->first,
    WProgenitor = tree->outgoingLines().rbegin()->first;
  if(abs(WProgenitor->id())!=ParticleID::Wplus) 
    swap(bProgenitor,WProgenitor);
  // Get the top quark
  ShowerProgenitorPtr topProgenitor = tree->incomingLines().begin()->first;
  // masses of the particles
  mt_ = topProgenitor->progenitor()->momentum().mass();
  w_  = WProgenitor  ->progenitor()->momentum().mass() / mt_;
  b_  = bProgenitor  ->progenitor()->momentum().mass() / mt_; // find rotation fgrom lab to frame with W along -z
  LorentzRotation eventFrame( topProgenitor->progenitor()->momentum().findBoostToCM() );
  Lorentz5Momentum pspectator = eventFrame*WProgenitor->progenitor()->momentum();
  eventFrame.rotateZ( -pspectator.phi() );
  eventFrame.rotateY( -pspectator.theta() - Constants::pi );



//   cerr << "testing " << *topProgenitor->progenitor() << "\n";
//   cerr << "testing " << *bProgenitor  ->progenitor() << "\n";
//   cerr << "testing " << *WProgenitor  ->progenitor() << "\n";


//   cerr << "testing " << eventFrame*(topProgenitor->progenitor()->momentum())/GeV << "\n";
//   cerr << "testing " << eventFrame*(bProgenitor  ->progenitor()->momentum())/GeV << "\n";
//   cerr << "testing " << eventFrame*(WProgenitor  ->progenitor()->momentum())/GeV << "\n";


  // invert it
  eventFrame.invert();
  // generate the hard emission
  vector<Lorentz5Momentum> momenta = hardMomenta();
  // if no emission return
  if(momenta.empty()) {
    ofstream file("empty.top", ios::app);
    file << "empty" << endl;
    file.close();
    topProgenitor->maximumpT(pTmin_);
    bProgenitor  ->maximumpT(pTmin_);
    return HardTreePtr();
  }
  // rotate momenta back to the lab
//   cerr << "testing size " << momenta.size() << "\n";
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    momenta[ix] *= eventFrame;
//     cerr << "new " 
// 	 << momenta[ix]/GeV << " " << momenta[ix].mass()/GeV << " "<<  momenta[ix].m()/GeV 
// 	 << "\n";
  }
  // get ParticleData objects
  tcPDPtr top    = topProgenitor->progenitor()->dataPtr();
  tcPDPtr bottom = bProgenitor  ->progenitor()->dataPtr();
  tcPDPtr Wboson = WProgenitor  ->progenitor()->dataPtr();
  tcPDPtr gluon  = getParticleData(ParticleID::g);
  // create new ShowerParticles
  ShowerParticlePtr emitter  (new_ptr(ShowerParticle(bottom,true )));
  ShowerParticlePtr spectator(new_ptr(ShowerParticle(Wboson,true )));
  ShowerParticlePtr gauge    (new_ptr(ShowerParticle(gluon ,true )));
  ShowerParticlePtr incoming (new_ptr(ShowerParticle(top   ,false)));
  ShowerParticlePtr parent   (new_ptr(ShowerParticle(bottom,true )));
  // set momenta
  emitter  ->set5Momentum(momenta[1]); 
  spectator->set5Momentum(momenta[2]);  
  gauge    ->set5Momentum(momenta[3]); 
  incoming ->set5Momentum(topProgenitor->progenitor()->momentum());  
  Lorentz5Momentum parentMomentum(momenta[1]+momenta[3]);
  parentMomentum.rescaleMass();
  parent->set5Momentum(parentMomentum);
  // Create the vectors of HardBranchings to create the HardTree:
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  // Incoming top quark
  spaceBranchings.push_back(new_ptr(HardBranching(incoming,SudakovPtr(),
						  HardBranchingPtr(),
						  HardBranching::Incoming)));
  // Outgoing particles from hard emission:
  HardBranchingPtr spectatorBranch(new_ptr(HardBranching(spectator,SudakovPtr(),
							 HardBranchingPtr(),
							 HardBranching::Outgoing)));
  HardBranchingPtr emitterBranch(new_ptr(HardBranching(parent,SudakovPtr(),
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(emitter,SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(gauge,SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  allBranchings.push_back(spaceBranchings[0]);
  allBranchings.push_back(emitterBranch);
  allBranchings.push_back(spectatorBranch);
  // Make the HardTree from the HardBranching vectors.
  HardTreePtr hardtree = new_ptr(HardTree(allBranchings,spaceBranchings,
					  ShowerInteraction::QCD));
  // Set the maximum pt for all other emissions
  topProgenitor->maximumpT(pT_);
  bProgenitor  ->maximumpT(pT_);
  // Connect the particles with the branchings in the HardTree
  hardtree->connect( topProgenitor->progenitor(), spaceBranchings[0] );
  hardtree->connect(   bProgenitor->progenitor(),   allBranchings[1] );
  hardtree->connect(   WProgenitor->progenitor(),   allBranchings[2] );
  // colour flow
  ColinePtr newline=new_ptr(ColourLine());
  for(set<HardBranchingPtr>::const_iterator cit=hardtree->branchings().begin();
      cit!=hardtree->branchings().end();++cit) {
    if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3)
      newline->addColoured((**cit).branchingParticle());
    else if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3bar)
      newline->addAntiColoured((**cit).branchingParticle());
  }

//   cerr << *hardtree << "\n";
//   exit(0);
  // return the tree
  return hardtree;
}


vector<Lorentz5Momentum>  SMTopPOWHEGDecayer::hardMomenta() {
  vector<Lorentz5Momentum> particleMomenta (4);
  double ymin = -10.;
  double ymax = 10.;
  double C = 4.;
  Energy2 lambda = sqr(mt_)* sqrt( 1. + pow(w_,4) + pow(b_,4) - 
			         2.*sqr(w_) - 2.*sqr(b_) - 2.*sqr(w_*b_));    

  //Calculate A
  double A = (ymax - ymin) * C * (coupling()->overestimateValue() / 
				 (2.*Constants::pi));
  Energy pTmax = mt_* (sqr(1.-w_) - sqr(b_)) / (2.*(1.-w_));

  while (pTmax >= pTmin_){  
    //Generate pT, y and phi values
    Energy pT = pTmax * pow(UseRandom::rnd() , (1./A));  
    if (pT < pTmin_) {
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
      InvEnergy2 meRatio = matrixElementRatio (particleMomenta);
      
      //Calculate jacobian 
      double J = (sqr(mt_) * particleMomenta[2].vect().mag2()) / 
	         (8. * pow(Constants::pi,3) * lambda * 
		 (particleMomenta[2].vect().mag()*(mt_-particleMomenta[3].e()) -
		  particleMomenta[2].e()*particleMomenta[3].z()));
      
      //Calculate weight
      weight[j] = sqr(pT) * meRatio * fabs(J) * coupling()->ratio(pT*pT) / C;
    }

    ofstream weights;
    if (weight[0] + weight[1] > 1.){
      weights.open("weights.top", ios::app);
      weights << weight[0]+weight[1] << endl;
    }

    //Accept point if weight > R
    if (weight[0] + weight[1] > UseRandom::rnd()) {
      if (weight[0] > (weight[0] + weight[1])*UseRandom::rnd()) {
	particleMomenta[2].setE((mt_/2.)*xw[0]);
	particleMomenta[2].setZ(-(mt_/2.)*sqrt(sqr(xw[0])-4.*sqr(w_)));
      }
      else {
	particleMomenta[2].setE((mt_/2.)*xw[1]);
	particleMomenta[2].setZ(-(mt_/2.)*sqrt(sqr(xw[1])-4.*sqr(w_)));
      }
      pT_ = pT;
      break;   
    }
    //If there's no splitting lower the pT
    pTmax = pT; 
  }
  return particleMomenta;
}


InvEnergy2 SMTopPOWHEGDecayer::matrixElementRatio(
			   vector<Lorentz5Momentum> particleMomenta) {
             
  Energy2 f = sqr(mt_) * (1. + pow(b_,4) - 2.*pow(w_,4) + 
			sqr(w_) + sqr(w_*b_) - 2.*sqr(b_));
  double Nc = standardModel()->Nc();
  double Cf = (sqr(Nc) - 1.) / (2.*Nc);  
  Energy2 B = (1./(2.*sqr(w_)))*f;  
  double Norm = sqr(Constants::pi)*Cf/2.;

  Energy2 PbPg = particleMomenta[1]*particleMomenta[3];
  Energy2 PtPg = particleMomenta[0]*particleMomenta[3];
  Energy2 PtPb = particleMomenta[0]*particleMomenta[1];

  double R = Norm * ( (-4.*f/sqr(w_)) * 
          ((sqr(mt_*b_)/sqr(PbPg)) + (sqr(mt_)/sqr(PtPg)) -2.*(PtPb/(PbPg*PtPg))) +
          (16. + 8.*sqr(1./w_) + 8.*sqr(b_/w_))* ((PtPg/PbPg) + (PbPg/PtPg)) -
	  (16./sqr(w_)) * (1. + sqr(b_)) ); 
   
  return R/B;
}


bool SMTopPOWHEGDecayer::calcMomenta(int j, Energy pT, double y, double phi,
				     double& xg, double& xw,
				     vector<Lorentz5Momentum>& particleMomenta){
  //Calculate xg
  xg = 2.*pT*cosh(y) / mt_;
  if (xg>(1. - sqr(b_ + w_)) || xg<0.) return false;

  //Calculate xw
  double xT = 2.*pT / mt_;
  double A = 4. - 4.*xg + sqr(xT);
  double B = 4.*(3.*xg - 2. + 2.*sqr(b_) - 2.*sqr(w_) - sqr(xg) - 
		 xg*sqr(b_) + xg*sqr(w_));
  double L = 1. + pow(w_,4) + pow(b_,4) - 2.*sqr(w_) - 2.*sqr(b_) - 2.*sqr(w_*b_);
  double det = 16.*(-L*sqr(xT) + pow(xT,4)*sqr(w_) + 
		    2.*xg*sqr(xT)*(1. - sqr(w_) - sqr(b_)) + 
		    L*sqr(xg) - sqr(xg*xT)*(1. + sqr(w_)) + 
		    2.*pow(xg,3)*(- 1. + sqr(w_) + sqr(b_)) + pow(xg,4));

  if (det<0.) return false;
  if (j==0) xw = (-B + sqrt(det))/(2.*A);
  if (j==1) xw = (-B - sqrt(det))/(2.*A);  
  if (xw>(1. + sqr(w_) - sqr(b_)) || xw<0.) return false;

  //Calculate xb
  double xb = 2. - xw - xg;     
  if (xb>(1. + sqr(b_) - sqr(w_)) || xb<0.) return false;       

  //Calculate xb_z  
  double xb_z;
  double epsilon_p =  -sqrt(sqr(xw) - 4.*sqr(w_)) + xT*sinh(y) +
                       sqrt(sqr(xb) - 4.*sqr(b_) - sqr(xT));
  double epsilon_m =  -sqrt(sqr(xw) - 4.*sqr(w_)) + xT*sinh(y) - 
                       sqrt(sqr(xb) - 4.*sqr(b_) - sqr(xT));

  if (fabs(epsilon_p) < 1.e-6){
    xb_z = sqrt(sqr(xb) - 4.*sqr(b_) - sqr(xT));
  }
  else if (fabs(epsilon_m) < 1.e-6){
    xb_z = -sqrt(sqr(xb) - 4.*sqr(b_) - sqr(xT));
  }
  else return false;

  //Check b is on shell
  if (fabs((sqr(xb) - sqr(xT) - sqr(xb_z) - 4.*sqr(b_)))>1.e-6) return false;

  //Calculate 4 momenta
  particleMomenta[0].setE(mt_);
  particleMomenta[0].setMass(mt_);

  particleMomenta[1].setE((mt_/2.)*xb);
  particleMomenta[1].setX(-pT*cos(phi));
  particleMomenta[1].setY(-pT*sin(phi));
  particleMomenta[1].setZ((mt_/2.)*xb_z);
  particleMomenta[1].setMass(mt_*b_);

  particleMomenta[2].setE((mt_/2.)*xw);
  particleMomenta[2].setZ(-(mt_/2.)*sqrt(sqr(xw) - 4.*sqr(w_)));
  particleMomenta[2].setMass(mt_*w_);

  particleMomenta[3].setE(pT*cosh(y));
  particleMomenta[3].setX(pT*cos(phi));
  particleMomenta[3].setY(pT*sin(phi));
  particleMomenta[3].setZ(pT*sinh(y));
 
  return true;
}


bool SMTopPOWHEGDecayer::psCheck(double xg, double xw) {
            
  //Check is point is in allowed region of phase space
  double xb_star = (1. - sqr(w_) + sqr(b_) - xg) / sqrt(1. - xg);
  double xg_star = xg / sqrt(1. - xg);

  if ((sqr(xb_star) - 4.*sqr(b_)) < 0) return false;
  double xw_max = (4. + 4.*sqr(w_) - sqr(xb_star + xg_star) + 
		   sqr(sqrt(sqr(xb_star) - 4.*sqr(b_)) + xg_star)) / 4.;
  double xw_min = (4. + 4.*sqr(w_) - sqr(xb_star + xg_star) + 
		   sqr(sqrt(sqr(xb_star) - 4.*sqr(b_)) - xg_star)) / 4.;

  if (xw < xw_min || xw > xw_max) return false;
  return true;
}

