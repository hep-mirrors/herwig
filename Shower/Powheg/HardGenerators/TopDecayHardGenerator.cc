// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TopDecayHardGenerator class.
//

#include "TopDecayHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"

using namespace Herwig;

TopDecayHardGenerator::TopDecayHardGenerator() 
  : pTmin_(1*GeV), overEstimate_(1000.), power_(1.0)
{}

IBPtr TopDecayHardGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr TopDecayHardGenerator::fullclone() const {
  return new_ptr(*this);
}

void TopDecayHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << alphaS_ << ounit(pTmin_,GeV) << overEstimate_ << power_
     << FFWVertex_ << FFGVertex_;
}

void TopDecayHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> alphaS_ >> iunit(pTmin_,GeV) >> overEstimate_ >> power_
     >> FFWVertex_ >> FFGVertex_;
}

ClassDescription<TopDecayHardGenerator> TopDecayHardGenerator::initTopDecayHardGenerator;
// Definition of the static class description member.

void TopDecayHardGenerator::Init() {

  static ClassDocumentation<TopDecayHardGenerator> documentation
    ("There is no documentation for the TopDecayHardGenerator class");

  static Reference<TopDecayHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &TopDecayHardGenerator::alphaS_, false, false, true, false, false);
  
  static Parameter<TopDecayHardGenerator, Energy> interfacePtMin
    ("minPt",
     "The pt cut on hardest emision generation",
     &TopDecayHardGenerator::pTmin_, GeV, 1.*GeV, 0*GeV, 100000.0*GeV,
     false, false, Interface::limited);

  static Parameter<TopDecayHardGenerator,double> interfaceOverEstimate
    ("OverEstimate",
     "The over estiamte for the veto algorithm",
     &TopDecayHardGenerator::overEstimate_, 1., 0., 10.0,
     false, false, Interface::limited);

  static Parameter<TopDecayHardGenerator,double> interfacePower
    ("Power",
     "The power for the sampling",
     &TopDecayHardGenerator::power_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

}

void TopDecayHardGenerator::doinit() {
  HardestEmissionGenerator::doinit();
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  // do the initialisation
  if(!hwsm) throw InitException() 
	      << "Must be the Herwig++ StandardModel class in "
	      << "TopDecayHardGenerator::doinit" << Exception::abortnow;
  FFWVertex_ = hwsm->vertexFFW();
  FFGVertex_ = hwsm->vertexFFG();
  // test of the phase-space volume
//   Energy Q  = getParticleData(ParticleID::t)->mass();
//   Energy mb = getParticleData(ParticleID::b)->mass();
//   Energy mw = getParticleData(ParticleID::Wplus)->mass();
//   double mu1 = mb/Q;
//   double mu2 = mw/Q;
//   //double mu1=0.;
//   //  double mu2=0.;
//   double wgt,wgtsum(0.),wgtsq(0.),wgtmax(1.);
//   ofstream output("dalitz.top");
//   output << "NEW FRAME\n";
//   output << "SET LIMITS X " << 0. << " " << 1.-sqr(mu1+mu2) << "\n";
//   output << "SET LIMITS Y " << 2*mu2 << " " << 1.+sqr(mu2)-sqr(mu1) << "\n"; 
//   unsigned long npoint(100000000);
//   for(unsigned long ix=0;ix<npoint;++ix) {
//     if(ix%50000==0) output << "PLOT\n";
//     double x1max = 1.+sqr(mu1)-sqr(mu2);
//     double x1min = 2.*mu1;
//     double x1 = x1min+UseRandom::rnd()*(x1max-x1min);
//     double lam = sqr(1.+sqr(mu1)-x1) + pow(mu2,4) - 2.*sqr(mu2)*(1.+sqr(mu1)-x1);
//     double x2min = 0.5/(1.+sqr(mu1)-x1)*((2.-x1)*(1+sqr(mu1)+sqr(mu2)-x1)
// 					 -sqrt((sqr(x1)-4.*sqr(mu1))*lam));
//     double x2max = 0.5/(1.+sqr(mu1)-x1)*((2.-x1)*(1+sqr(mu1)+sqr(mu2)-x1)
// 					 +sqrt((sqr(x1)-4.*sqr(mu1))*lam));
//     double x2 = x2min+UseRandom::rnd()*(x2max-x2min);
//     wgt = (x1max-x1min)*(x2max-x2min);
//     wgtsum += wgt;
//     wgtsq += sqr(wgt);
//     if(wgt<0) cerr << "negative weight \n";
//     if(wgt>wgtmax*UseRandom::rnd()&&ix<200000) output << 2-x1-x2 << "\t" << x2  << "\n";
//   }
//   output << "PLOT\n";
//   wgtsum = wgtsum/float(npoint);
//   wgtsq  = max(wgtsq/float(npoint)-sqr(wgtsum),0.);
//   cerr << "testing int A " << wgtsum << " +/- " << sqrt(wgtsq/float(npoint)) << "\n"; 
//   cerr << "max weight " << wgtmax << "\n";
//   wgtsum = wgtsq = wgtmax = 0.;
//   output << "NEW FRAME\n";
//   output << "SET LIMITS X " << 0. << " " << 1.-sqr(mu1+mu2) << "\n";
//   output << "SET LIMITS Y " << 2*mu2 << " " << 1.+sqr(mu2)-sqr(mu1) << "\n"; 
//   wgtmax = 150.;
//   npoint *=10;
//   for(unsigned long ix=0;ix<npoint;++ix) {
//     if(ix%50000==0) output << "PLOT\n";
//     Energy ptmax = 0.5*Q;
//     Energy pt = UseRandom::rnd()*ptmax;
//     double ymax = 0.5*Q/pt*(1.-sqr(mu1+mu2));
//     if(ymax<1.) continue;
//     ymax = acosh(ymax);
//     double y   = ymax*(2.*UseRandom::rnd()-1);
//     double x3 = 2.*pt/Q*cosh(y),x1,x2;
//     // check x3 within limits
//     if(x3<0||x3>1.-sqr(mu1+mu2)) continue;
//     // possible values of x2
//     double r = 4.*sqr(pt/Q)*sqr(sinh(y));
//     double A = sqr(2.-x3)-r;
//     double C = -2.*(1.-x3-sqr(mu1)+sqr(mu2));
//     double B = 2.*(2.-x3)*C;
//     C = sqr(C)+4.*r*sqr(mu2);
//     double det = sqr(B)-4.*A*C;
//     if(det<0.) continue;
//     det = sqrt(det);
//     double test[2] = { 0.5/A*(-B+det) , 0.5/A*(-B-det) };
//     // pick the one which gives the right sign of y
//     for(unsigned int iy=0;iy<2;++iy) {
//       x2 = test[iy];
//       double ratio = 0.5*Q/sqrt(sqr(x2)-4.*sqr(mu2))*
// 	(x2*(2.-x3)-2.*(1.-x3-sqr(mu1)+sqr(mu2)))/pt/sinh(y);
//       if(ratio>0.) break;
//     }
//     // check within the limits
//     double lam = sqrt(sqr(1.-x3)+pow(mu1,4)+pow(mu2,4)-2.*sqr(mu1*mu2)
// 		      -2.*(1.-x3)*(sqr(mu2)+sqr(mu1)));
//     double lim[2] = {0.5/(1.-x3)*((2.-x3)*(1.-x3-sqr(mu1)+sqr(mu2))-x3*lam),
// 		     0.5/(1.-x3)*((2.-x3)*(1.-x3-sqr(mu1)+sqr(mu2))+x3*lam)};
//     if(x2<lim[0]||x2>lim[1]) continue;
//     // finally x1
//     x1 = 2. - x2 - x3;
//     Energy2 J =-0.5*sqr(Q)*(x2*sqr(mu2)+2*x1*sqr(mu2)-sqr(x2)+x2-x2*x1+x2*sqr(mu1))/
//       pow(sqrt(sqr(x2)-4*sqr(mu2)),3);
//     wgt = pt*ptmax*2.*ymax/abs(J);
//     if(wgt>wgtmax*UseRandom::rnd()&&ix<1000000) output << 2-x1-x2 << "\t" << x2  << "\n";
//     wgtsum += wgt;
//     wgtsq  += sqr(wgt);
//     if(wgt<0) cerr << "negative weight \n";
//   }
//   output << "PLOT\n";
//   output.close();
//   wgtsum = wgtsum/float(npoint);
//   wgtsq  = max(wgtsq/float(npoint)-sqr(wgtsum),0.);
//   cerr << "testing int B " << wgtsum << " +/- " << sqrt(wgtsq/float(npoint)) << "\n";
//   cerr << "max weight " << wgtmax << "\n";
}

HardTreePtr TopDecayHardGenerator::generateHardest(ShowerTreePtr tree) {
  vector<ShowerProgenitorPtr> particlesToShower;
  // if no emission limit emission in shower
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit = tree->incomingLines().begin(); 
      cit != tree->incomingLines().end(); ++cit ) {
    particlesToShower.push_back(cit->first);
    cit->first->maximumpT(pTmin_);
  }
  for( map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	 cit= tree->outgoingLines().begin();
       cit != tree->outgoingLines().end(); ++cit ) {
    particlesToShower.push_back(cit->first);
    cit->first->maximumpT(pTmin_);
  }
  tcPDPtr gluon = getParticleData(ParticleID::g);
  // extract the top mass
  Energy mt = tree->incomingLines().begin()->first->progenitor()->mass();
  pair<tcPDPtr,Lorentz5Momentum> top = 
    make_pair(tree->incomingLines().begin()->first->progenitor()->dataPtr(),
	      tree->incomingLines().begin()->first->progenitor()->momentum());
  // extract the spin density matrix for the top
  tSpinPtr sp = tree->incomingLines().begin()->first->copy()->spinInfo();
  tFermionSpinPtr inspin = !sp ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(sp);
  RhoDMatrix rho = inspin ? inspin->rhoMatrix() : RhoDMatrix(PDT::Spin1Half);
  // bottom and W mass
  Energy mb,mw;
  pair<tcPDPtr,Lorentz5Momentum> bottom,Wboson;
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator it=
	tree->outgoingLines().begin(); it != tree->outgoingLines().end(); ++it) {
    long id = abs(it->first->progenitor()->id());
    if(id==ParticleID::b)          {
      mb = it->first->progenitor()->mass();
      bottom = make_pair(it->first->progenitor()->dataPtr(),
			 it->first->progenitor()->momentum());
    }
    else if(id==ParticleID::Wplus) {
      mw = it->first->progenitor()->mass();
      Wboson = make_pair(it->first->progenitor()->dataPtr(),
			 it->first->progenitor()->momentum());
    }
  }
  double mub = mb/mt, muw = mw/mt;
  double lam2 = sqrt((1.-sqr(muw+mub))*(1.-sqr(muw-mub)));
  // W decay products
  vector<pair<tcPDPtr,Lorentz5Momentum> > Wdec;
  LorentzRotation Wboost(-tree->treelinks().begin()->first->
			 incomingLines().begin()->second->momentum().boostVector());
  Wboost.boost(Wboson.second.boostVector());
  if(tree->treelinks().size()==1) {
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  it=tree->treelinks().begin()->first->outgoingLines().begin();
	it!=tree->treelinks().begin()->first->outgoingLines().end();++it) {
      Lorentz5Momentum pnew = Wboost*it->second->momentum();
      Wdec.push_back(make_pair(it->second->dataPtr(),pnew));
    }
    if(Wdec[0].first->id()<0) swap(Wdec[0],Wdec[1]);
  }
  // compute the leading-order matrix element
  double lome = loME(rho,top,bottom,Wboson,Wdec);
  // rotation for new momenta
  LorentzRotation rot(-top.second.boostVector());
  Axis axis((rot*bottom.second).vect().unit());
  if(axis.perp2()>0.) {
    double sinth(-sqrt(sqr(axis.x())+sqr(axis.y())));
    rot.rotate(acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  }
  else if(axis.z()<0 ) {
    rot.rotateX(Constants::pi);
  }
  rot.invert();
  // generate the pT and rapidity of the emission
  Energy pT(mt);
  double ymax = acosh(0.5*mt/pTmin_*(1.-sqr(mub+muw)));
  double a = 2.*overEstimate_*ymax*alphaS_->overestimateValue()/Constants::twopi;
  double y,phi;
  double x1,x2,x3;
  Lorentz5Momentum pw,pg,pb;
  do {
    // generate kinematic variables
    if(power_== 1. ) {
      pT *= pow(UseRandom::rnd(),1./a);
    }
    else {
      pT = mt*pow((1.-power_)/a*log(UseRandom::rnd())+pow(pT/mt,1.-power_),1./(1.-power_));
    }
    y   = ymax*(2.*UseRandom::rnd()-1);
    phi = UseRandom::rnd()*Constants::twopi;
    // kinematic variables
    x3 = 2.*pT/mt*cosh(y);
    // check x3 within limits
    if(x3<0||x3>1.-sqr(mub+muw)) continue;
    // possible values of x2
    double r = 4.*sqr(pT/mt)*sqr(sinh(y));
    double A = sqr(2.-x3)-r;
    double C = -2.*(1.-x3-sqr(mub)+sqr(muw));
    double B = 2.*(2.-x3)*C;
    C = sqr(C)+4.*r*sqr(muw);
    double det = sqr(B)-4.*A*C;
    if(det<0.) continue;
    det = sqrt(det);
    // two solns
    if(abs(det/B)>1e-8) {
      // pick the one which gives the right sign of y
      double test[2] = { 0.5/A*(-B+det) , 0.5/A*(-B-det) };
      for(unsigned int ix=0;ix<2;++ix) {
	x2 = test[ix];
	double ratio = 0.5*mt/sqrt(sqr(x2)-4.*sqr(muw))*
	  (x2*(2.-x3)-2.*(1.-x3-sqr(mub)+sqr(muw)))/pT/sinh(y);
	if(ratio>0.) break;
      }
    }
    // only one soln
    else {
      x2 =  -0.5/A*B;
    }
    // check within the limits
    double lam = sqrt(sqr(1.-x3)+pow(mub,4)+pow(muw,4)-2.*sqr(mub*muw)
		      -2.*(1.-x3)*(sqr(muw)+sqr(mub)));
    double lim[2] = {0.5/(1.-x3)*((2.-x3)*(1.-x3-sqr(mub)+sqr(muw))-x3*lam),
		     0.5/(1.-x3)*((2.-x3)*(1.-x3-sqr(mub)+sqr(muw))+x3*lam)};
    if(x2<lim[0]||x2>lim[1]) continue;
    // finally x1
    x1 = 2. - x2 - x3;
    // compute the new momenta
    pw = Lorentz5Momentum(ZERO,ZERO,-0.5*mt*sqrt(sqr(x2)-4.*sqr(muw)),0.5*mt*x2);
    pg = Lorentz5Momentum(pT*cos(phi),pT*sin(phi),pT*sinh(y),pT*cosh(y));
    pb = Lorentz5Momentum(ZERO,ZERO,ZERO,mt)-pw-pg;
    pw.rescaleMass();
    pb.rescaleMass();
    pg.rescaleMass();
    pw *= rot;
    pb *= rot;
    pg *= rot;
    // boost for the leptons if needed
    vector<pair<tcPDPtr,Lorentz5Momentum> > newdec=Wdec;
    if(!Wdec.empty()) {
      double sx = 0.5/sqr(mt)/sqr(mw)*
	((sqr(mt)+sqr(mw)-sqr(mb))*sqrt(sqr(top.second*pw)-sqr(mw*mt))
	 -(top.second*pw)*sqrt((sqr(mt)-sqr(mb-mw))*(sqr(mt)-sqr(mb+mw))));
      double cx = sqrt(1+sqr(sx));
      Energy2 D = sqrt(sqr(top.second*pw)-sqr(mw*mt));
      for(unsigned int ix=0;ix<newdec.size();++ix) {
	newdec[ix].second +=
	  -sx/D*(top.second*(pw*newdec[ix].second)-pw*(top.second*newdec[ix].second))
	  +(cx-1)/sqr(D)*((top.second*pw)*(top.second*(pw*newdec[ix].second)
					   +pw*(top.second*newdec[ix].second))
			  -sqr(mw)*top.second*(top.second*newdec[ix].second)
			  -sqr(mt)*pw*(pw*newdec[ix].second));
      }
    }
    // compute the ratio of the real emission to leading order matrix element
    pair<tcPDPtr,Lorentz5Momentum> newW = make_pair(Wboson.first,pw);
    pair<tcPDPtr,Lorentz5Momentum> newB = make_pair(bottom.first,pb);
    pair<tcPDPtr,Lorentz5Momentum> newG = make_pair(gluon       ,pg);
    // compute the rejection weight
    // matrix element weight
    double meratio = realME(rho,top,newB,newG,newW,newdec)/lome;
    // jacobian from x1,x2 to pT,y
    Energy2 J = abs(0.5*sqr(mt)*(x2*sqr(muw)+2*x1*sqr(muw)-sqr(x2)+x2-x2*x1+x2*sqr(mub))/
		    pow(sqr(x2)-4*sqr(muw),1.5));
    // total weight
    double wgt = 0.5/overEstimate_*alphaS_->ratio(sqr(pT))*mt*pT/J/lam2*meratio; 
    if(power_==1) wgt *= pT/mt;
    else          wgt *= pow(pT/mt,power_);
    // check for wgt>1
    if(wgt>1.) cerr << "testing weight problem " << wgt << "\n";
    // rejection step
    if(wgt>UseRandom::rnd()) break;
  }
  while (pT>pTmin_);
  // return if no emission
  if(pT<pTmin_) return HardTreePtr();
  // sudakov for emission from the bottom
  SudakovPtr sudakov;
  BranchingList branchings = 
    evolver()->splittingGenerator()->finalStateBranchings();
  long index = abs(bottom.first->id());
  for(BranchingList::const_iterator cit = branchings.lower_bound(index);
      cit != branchings.upper_bound(index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0] == index && ids[1] == index && ids[2] == gluon->id() ) {
      sudakov = cit->second.first;
      break; 	    
    }
  }
  ShowerParticlePtr newB  ( new_ptr( ShowerParticle( bottom.first, true  ) ) );
  ShowerParticlePtr newW  ( new_ptr( ShowerParticle( Wboson.first, true  ) ) );
  ShowerParticlePtr newG  ( new_ptr( ShowerParticle( gluon       , true  ) ) );
  ShowerParticlePtr newT  ( new_ptr( ShowerParticle( top.first   , false ) ) );
  ShowerParticlePtr parent( new_ptr( ShowerParticle( bottom.first, true  ) ) );
  newB->set5Momentum(pb);
  newW->set5Momentum(pw);
  newG->set5Momentum(pg);
  newT->set5Momentum(top.second);
  Lorentz5Momentum poff(pb+pg);
  poff.rescaleMass();
  parent->set5Momentum(poff);
  // Create the vectors of HardBranchings to create the HardTree:
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  // incoming top
  spaceBranchings.push_back(new_ptr(HardBranching(newT,SudakovPtr(),
						  HardBranchingPtr(),
						  HardBranching::Incoming)));
  // W
  HardBranchingPtr wBranch(new_ptr (HardBranching(newW,SudakovPtr(), 
						  HardBranchingPtr(),
						  HardBranching::Outgoing)));
  // bottom stuff
  HardBranchingPtr emitterBranch(new_ptr(HardBranching(parent,sudakov,
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(newB,SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(newG,SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  // add branchings to the tree
  allBranchings.push_back( emitterBranch );
  allBranchings.push_back( wBranch );
  allBranchings.push_back( spaceBranchings.back() );
  // Make the HardTree from the HardBranching vectors.
  HardTreePtr nasontree = new_ptr( HardTree( allBranchings, spaceBranchings ,
					     ShowerInteraction::QCD ) );
  // Connect the particles with the branchings in the HardTree
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    long id = particlesToShower[ix]->progenitor()->id();
    if(abs(id)==ParticleID::t) 
      nasontree->connect( particlesToShower[ix]->progenitor(), allBranchings[2] );
    else if(abs(id)==ParticleID::b) 
      nasontree->connect( particlesToShower[ix]->progenitor(), allBranchings[0] );
    else if(abs(id)==ParticleID::Wplus) 
      nasontree->connect( particlesToShower[ix]->progenitor(), allBranchings[1] );
  }
  // create the colour lines
  ColinePtr newLine  = new_ptr( ColourLine() );
  if(top.first->id()>0) {
    newLine->addColoured(newT);
    newLine->addColoured(parent);
  }
  else {
    newLine->addAntiColoured(newT);
    newLine->addAntiColoured(parent);
  }
  // Calculate the shower variables:
  evolver()->showerModel()->kinematicsReconstructor()->
    deconstructDecayJets(nasontree,evolver(),ShowerInteraction::QCD);
  // Return the HardTree
  return nasontree;
}

bool TopDecayHardGenerator::canHandle(ShowerTreePtr tree) {
  // require one incoming top
  if(tree->incomingLines().size()!=1) return false; 
  if(abs(tree->incomingLines().begin()->first->id()) != ParticleID::t) return false;
  if(tree->outgoingLines().size()!=2) return false;
  bool found[2]={false,false};
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator it=
	tree->outgoingLines().begin(); it != tree->outgoingLines().end(); ++it) {
    long id = abs(it->first->progenitor()->id());
    if(id==ParticleID::b) found[0] = true;
    else if(id==ParticleID::Wplus) found[1] = true;
  }
  return found[0] && found[1];
}

double TopDecayHardGenerator::loME(RhoDMatrix rho,
				   pair<tcPDPtr,Lorentz5Momentum> top, 
				   pair<tcPDPtr,Lorentz5Momentum> bottom, 
				   pair<tcPDPtr,Lorentz5Momentum> Wboson,
				   vector<pair<tcPDPtr,Lorentz5Momentum> > leptons) {
  Energy2 scale = top.second.mass2();
  // wavefunctions for the top/bottom, don't depend on 2/3 body
  SpinorWaveFunction    wave1;
  SpinorBarWaveFunction wave2;
  if(top.first->id()>0) {
    wave1 = SpinorWaveFunction   (   top.second,   top.first,incoming);
    wave2 = SpinorBarWaveFunction(bottom.second,bottom.first,outgoing);
  }
  else {
    wave2 = SpinorBarWaveFunction(   top.second,   top.first,incoming);
    wave1 = SpinorWaveFunction   (bottom.second,bottom.first,outgoing);
  }
  // wavefunctions for the W
  vector<VectorWaveFunction> wave3;
  // two-body matrix element
  if(leptons.empty()) {
    VectorWaveFunction wwave(Wboson.second,Wboson.first,outgoing);
    for(unsigned int whel=0;whel<3;++whel) {
      wwave.reset(whel);
      wave3.push_back(wwave);
    }
  }
  // three-body matrix element
  else {
    SpinorBarWaveFunction fwave(leptons[0].second,leptons[0].first,outgoing);
    SpinorWaveFunction    awave(leptons[1].second,leptons[1].first,outgoing);
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      fwave.reset(ihel1);
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	awave.reset(ihel2);	  
	wave3.push_back(FFWVertex_->evaluate(scale,1,Wboson.first,awave,fwave));
      }
    }
  }
  Complex me[2][2][4];
  // compute the matrix element squared
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    wave1.reset(ihel1);
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      wave2.reset(ihel2);
      for(unsigned int ihel3=0;ihel3<wave3.size();++ihel3) {
	if(top.first->id()>0) {
	  me[ihel1][ihel2][ihel3] = FFWVertex_->evaluate(scale,wave1,wave2,wave3[ihel3]);
	}
	else {
	  me[ihel2][ihel1][ihel3] = FFWVertex_->evaluate(scale,wave1,wave2,wave3[ihel3]);
	}
      }
    }
  }
  Complex output(0.);
  for(unsigned int inhel1=0;inhel1<2;++inhel1) {
    for(unsigned int inhel2=0;inhel2<2;++inhel2) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	for(unsigned int ihel3=0;ihel3<wave3.size();++ihel3) {
	  output += rho(inhel1,inhel2)*me[inhel1][ihel2][ihel3]*
	    conj(me[inhel2][ihel2][ihel3]);
	}
      }
    }
  }
  // scale factor for 2 body
  if(leptons.empty()) {
    output /= scale*UnitRemoval::InvE2;
  }
  return output.real();
}

double TopDecayHardGenerator::realME(RhoDMatrix rho,
				     pair<tcPDPtr,Lorentz5Momentum> top, 
				     pair<tcPDPtr,Lorentz5Momentum> bottom, 
				     pair<tcPDPtr,Lorentz5Momentum> gluon ,
				     pair<tcPDPtr,Lorentz5Momentum> Wboson,
				     vector<pair<tcPDPtr,Lorentz5Momentum> > leptons) {
  Energy2 scale = top.second.mass2();
  // wavefunctions for the top/bottom, don't depend on 2/3 body
  SpinorWaveFunction    wave1;
  SpinorBarWaveFunction wave2;
  if(top.first->id()>0) {
    wave1 = SpinorWaveFunction   (   top.second,   top.first,incoming);
    wave2 = SpinorBarWaveFunction(bottom.second,bottom.first,outgoing);
  }
  else {
    wave2 = SpinorBarWaveFunction(   top.second,   top.first,incoming);
    wave1 = SpinorWaveFunction   (bottom.second,bottom.first,outgoing);
  }
  // wavefunctions for the gluon
  VectorWaveFunction wave4(gluon.second,gluon.first,outgoing);
  // wavefunctions for the W
  vector<VectorWaveFunction> wave3;
  // two-body matrix element
  if(leptons.empty()) {
    VectorWaveFunction wwave(Wboson.second,Wboson.first,outgoing);
    for(unsigned int whel=0;whel<3;++whel) {
      wwave.reset(whel);
      wave3.push_back(wwave);
    }
  }
  // three-body matrix element
  else {
    SpinorBarWaveFunction fwave(leptons[0].second,leptons[0].first,outgoing);
    SpinorWaveFunction    awave(leptons[1].second,leptons[1].first,outgoing);
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      fwave.reset(ihel1);
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	awave.reset(ihel2);	  
	wave3.push_back(FFWVertex_->evaluate(scale,1,Wboson.first,awave,fwave));
      }
    }
  }
  Complex me[2][2][2][4];
  // compute the matrix element squared
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    wave1.reset(ihel1);
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      wave2.reset(ihel2);
      for(unsigned int ig=0;ig<2;++ig) {
	wave4.reset(2*ig);
	// intermediates
	SpinorWaveFunction    inter1 = 
	  FFGVertex_->evaluate(scale,3,wave1.particle(),wave1,wave4);
	SpinorBarWaveFunction inter2 = 
	  FFGVertex_->evaluate(scale,3,wave2.particle(),wave2,wave4);
	for(unsigned int ihel3=0;ihel3<wave3.size();++ihel3) {
	  //  first diagram
	  Complex diag1 = FFWVertex_->evaluate(scale,inter1, wave2,wave3[ihel3]);
	  // second diagram
	  Complex diag2 = FFWVertex_->evaluate(scale, wave1,inter2,wave3[ihel3]);
	  if(top.first->id()>0) {
	    me[ihel1][ihel2][ig][ihel3] = diag1+diag2;
	  }
	  else {
	    me[ihel2][ihel1][ig][ihel3] = diag1+diag2;
	  }
	}
      }
    }
  }
  // sum over spins
  Complex output(0.);
  for(unsigned int inhel1=0;inhel1<2;++inhel1) {
    for(unsigned int inhel2=0;inhel2<2;++inhel2) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	for(unsigned int ig=0;ig<2;++ig) {
	  for(unsigned int ihel3=0;ihel3<wave3.size();++ihel3) {
	    output += rho(inhel1,inhel2)*me[inhel1][ihel2][ig][ihel3]*
	      conj(me[inhel2][ihel2][ig][ihel3]);
	  }
	}
      }
    }
  }
  // scale factor for 3 body
  if(!leptons.empty()) {
    output *= scale*UnitRemoval::InvE2;
  }
  // remove the strong coupling
  output /= generator()->standardModel()->alphaS(scale)*4.*Constants::pi;
  // colour factor
  output *= 4./3.;
  // return the answer
  return output.real();
}
