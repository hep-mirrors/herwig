// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISHardGenerator class.
//

#include "DISHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

DISHardGenerator::DISHardGenerator() : comptonWeight_(50.), BGFWeight_(150.), 
				       pTmin_(1.*GeV) 
{}

IBPtr DISHardGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr DISHardGenerator::fullclone() const {
  return new_ptr(*this);
}

void DISHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << alphaS_ << sinW_ << cosW_ << ounit(mz2_,GeV2) << ounit(pTmin_,GeV)
     << comptonWeight_ << BGFWeight_ << gluon_;
}

void DISHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> alphaS_ >> sinW_ >> cosW_ >> iunit(mz2_,GeV2) >> iunit(pTmin_,GeV)
     >> comptonWeight_ >> BGFWeight_ >> gluon_;
}

void DISHardGenerator::doinit() {
  HardestEmissionGenerator::doinit();
  // electroweak parameters
  sinW_ = generator()->standardModel()->sin2ThetaW();
  cosW_ = sqrt(1.-sinW_);
  sinW_ = sqrt(sinW_);
  mz2_ = sqr(getParticleData(ParticleID::Z0)->mass());
  gluon_ = getParticleData(ParticleID::g);
}

ClassDescription<DISHardGenerator> DISHardGenerator::initDISHardGenerator;
// Definition of the static class description member.

void DISHardGenerator::Init() {

  static ClassDocumentation<DISHardGenerator> documentation
    ("There is no documentation for the DISHardGenerator class");

  static Reference<DISHardGenerator,ShowerAlpha> interfaceShowerAlphaQCD
    ("ShowerAlphaQCD",
     "The object calculating the strong coupling constant",
     &DISHardGenerator::alphaS_, false, false, true, false, false);

  static Parameter<DISHardGenerator,Energy> interfacepTMin
    ("pTMin",
     "The minimum pT",
     &DISHardGenerator::pTmin_, GeV, 1.*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DISHardGenerator,double> interfaceComptonWeight
    ("ComptonWeight",
     "Weight for the overestimate ofthe compton channel",
     &DISHardGenerator::comptonWeight_, 50.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<DISHardGenerator,double> interfaceBGFWeight
    ("BGFWeight",
     "Weight for the overestimate of the BGF channel",
     &DISHardGenerator::BGFWeight_, 100.0, 0.0, 1000.0,
     false, false, Interface::limited);

}

bool DISHardGenerator::canHandle(ShowerTreePtr tree) {
  // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // two outgoing particles
  if(tree->outgoingLines().size()!=2) return false;
  // check incoming quark and lepton
  bool quark(false),lepton(false);
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    quark  |=  QuarkMatcher::Check(cit->first->progenitor()->data());
    lepton |= LeptonMatcher::Check(cit->first->progenitor()->data());
  }
  if(!quark||!lepton) return false;
  // check outgoing quark and lepton
  quark = false;
  lepton = false;
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    quark  |=  QuarkMatcher::Check(cjt->first->progenitor()->data());
    lepton |= LeptonMatcher::Check(cjt->first->progenitor()->data());
  }
  if(!quark||!lepton) return false;
  // can handle it
  return true;
}

HardTreePtr DISHardGenerator::generateHardest(ShowerTreePtr tree) {
  ShowerParticlePtr quark[2],lepton[2];
  PPtr hadron;
  // incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data())) {
      hadron = cit->first->original()->parents()[0];
      quark [0] = cit->first->progenitor();
      beam_ = cit->first->beam();
    }
    else if(LeptonMatcher::Check(cit->first->progenitor()->data())) {
      lepton[0] = cit->first->progenitor();
    }
  }
  pdf_=beam_->pdf();
  assert(beam_&&pdf_&&quark[0]&&lepton[0]);
  // outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data()))
      quark [1] = cit->first->progenitor();
    else if(LeptonMatcher::Check(cit->first->progenitor()->data()))
      lepton[1] = cit->first->progenitor();
  }
  assert(quark[1]&&lepton[1]);
  // Particle data objects
  for(unsigned int ix=0;ix<2;++ix) partons_[ix] = quark[ix]->dataPtr();
  // extract the born variables
  q_ =lepton[0]->momentum()-lepton[1]->momentum();
  q2_ = -q_.m2();
  xB_ = quark[0]->x();
  double  yB = 
    (                   q_*quark[0]->momentum())/
    (lepton[0]->momentum()*quark[0]->momentum()); 
  l_ = 2./yB-1.;
  // construct lorentz transform from lab to breit frame
  Lorentz5Momentum phadron =  hadron->momentum();
  phadron.setMass(0.*GeV);
  phadron.rescaleRho();
  Lorentz5Momentum pcmf = phadron+0.5/xB_*q_;
  pcmf.rescaleMass();
  rot_ = LorentzRotation(-pcmf.boostVector());
  Lorentz5Momentum pbeam = rot_*phadron;
  Axis axis(pbeam.vect().unit());
  double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
  rot_.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  Lorentz5Momentum pl    = rot_*lepton[0]->momentum();
  rot_.rotateZ(-atan2(pl.y(),pl.x()));
  // momenta of the particles
  pl_[0]=rot_*lepton[0]->momentum();
  pl_[1]=rot_*lepton[1]->momentum();
  pq_[0]=rot_* quark[0]->momentum();
  pq_[1]=rot_* quark[1]->momentum();
  q_ *= rot_;
  // coefficient for the matrix elements
  acoeff_ = A(lepton[0]->dataPtr(),lepton[1]->dataPtr(),
	      quark [0]->dataPtr(),quark [1]->dataPtr());
  // generate a compton point
  generateCompton();
  generateBGF();
  // no valid emission, return
  if(pTCompton_<ZERO&&pTBGF_<ZERO) return HardTreePtr();
  // type of emission, pick highest pT
  bool isCompton=pTCompton_>pTBGF_;
  // find the sudakov for the branching
  SudakovPtr sudakov;
  // ISR
  if(ComptonISFS_||!isCompton) {
    BranchingList branchings=evolver()->splittingGenerator()->initialStateBranchings();
    long index = abs(partons_[0]->id());
    IdList br(3);
    if(isCompton) {
      br[0] = index;
      br[1] = index;
      br[2] = ParticleID::g;
    }
    else {
      br[0] = ParticleID::g;
      br[1] =  abs(partons_[0]->id());
      br[2] = -abs(partons_[0]->id());
    }
    for(BranchingList::const_iterator cit = branchings.lower_bound(index); 
	cit != branchings.upper_bound(index); ++cit ) {
      IdList ids = cit->second.second;
      if(ids[0]==br[0]&&ids[1]==br[1]&&ids[2]==br[2]) {
	sudakov=cit->second.first;
	break;
      }
    }
  }
  // FSR
  else {
    BranchingList branchings = 
      evolver()->splittingGenerator()->finalStateBranchings();
    long index=abs(partons_[1]->id());
    for(BranchingList::const_iterator cit = branchings.lower_bound(index);
	cit != branchings.upper_bound(index); ++cit ) {
      IdList ids = cit->second.second;
      if(ids[0]==index&&ids[1]==index&&ids[2]==ParticleID::g) {
	sudakov = cit->second.first;
	break; 	    
      }
    }
  }
  if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
				 << "DISHardGenerator::generateHardest()" 
				 << Exception::runerror;
  // add the leptons
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  spaceBranchings.push_back(new_ptr(HardBranching(lepton[0],SudakovPtr(),
						  HardBranchingPtr(),
						  HardBranching::Incoming)));
  allBranchings.push_back(spaceBranchings.back());
  allBranchings.push_back(new_ptr(HardBranching(lepton[1],SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  // compton hardest
  if(isCompton) {
    rot_.invert();
    for(unsigned int ix=0;ix<ComptonMomenta_.size();++ix) {
      ComptonMomenta_[ix].transform(rot_);
    }
    ShowerParticlePtr newqout (new_ptr(ShowerParticle(partons_[1],true)));
    newqout->set5Momentum(ComptonMomenta_[1]);
    ShowerParticlePtr newg(new_ptr(ShowerParticle(gluon_,true)));
    newg->set5Momentum(ComptonMomenta_[2]);
    ShowerParticlePtr newqin   (new_ptr(ShowerParticle(partons_[0],false )));
    newqin->set5Momentum(ComptonMomenta_[0]);
    if(ComptonISFS_) {
      ShowerParticlePtr newspace(new_ptr(ShowerParticle(partons_[0],false)));
      newspace->set5Momentum(ComptonMomenta_[0]-ComptonMomenta_[2]);
      HardBranchingPtr spaceBranch(new_ptr(HardBranching(newqin,sudakov,
							 HardBranchingPtr(),
							 HardBranching::Incoming)));
      HardBranchingPtr offBranch(new_ptr(HardBranching(newspace,SudakovPtr(),
						       spaceBranch,
						       HardBranching::Incoming)));
      spaceBranch->addChild(offBranch);
      HardBranchingPtr g(new_ptr(HardBranching(newg,SudakovPtr(),spaceBranch,
					       HardBranching::Outgoing)));
      spaceBranch->addChild(g);
      HardBranchingPtr outBranch(new_ptr(HardBranching(newqout,SudakovPtr(),
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
      spaceBranchings.push_back(spaceBranch);
      allBranchings.push_back(offBranch);
      allBranchings.push_back(outBranch);
      ColinePtr newin(new_ptr(ColourLine())),newout(new_ptr(ColourLine()));
      newin->addColoured(newqin,partons_[0]->id()<0);
      newin->addColoured(newg  ,partons_[0]->id()<0);
      newout->addColoured(newspace,partons_[0]->id()<0);
      newout->addColoured(newqout,partons_[1]->id()<0);
      newout->addColoured(newg  ,partons_[1]->id()>0);
    }
    else {
      ShowerParticlePtr newtime(new_ptr(ShowerParticle(partons_[1],true)));
      newtime->set5Momentum(ComptonMomenta_[1]+ComptonMomenta_[2]);
      HardBranchingPtr spaceBranch(new_ptr(HardBranching(newqin,SudakovPtr(),
							 HardBranchingPtr(),
							 HardBranching::Incoming)));
      HardBranchingPtr offBranch(new_ptr(HardBranching(newtime,sudakov,
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
      HardBranchingPtr g(new_ptr(HardBranching(newg,SudakovPtr(),offBranch,
					       HardBranching::Outgoing)));
      HardBranchingPtr outBranch(new_ptr(HardBranching(newqout,SudakovPtr(),offBranch,
						       HardBranching::Outgoing)));
      offBranch->addChild(outBranch);
      offBranch->addChild(g);
      spaceBranchings.push_back(spaceBranch);
      allBranchings.push_back(spaceBranch);
      allBranchings.push_back(offBranch);
      ColinePtr newin(new_ptr(ColourLine())),newout(new_ptr(ColourLine()));
      newin->addColoured(newqin,partons_[0]->id()<0);
      newin->addColoured(newtime,partons_[0]->id()<0);
      newin->addColoured(newg  ,partons_[0]->id()<0);
      newout->addColoured(newqout,partons_[1]->id()<0);
      newout->addColoured(newg,partons_[1]->id()>0);
    }
  }
  // BGF hardest
  else {
    rot_.invert();
    for(unsigned int ix=0;ix<BGFMomenta_.size();++ix) {
      BGFMomenta_[ix].transform(rot_);
    }
    ShowerParticlePtr newq   (new_ptr(ShowerParticle(partons_[1],true)));
    newq->set5Momentum(BGFMomenta_[1]);
    ShowerParticlePtr newqbar(new_ptr(ShowerParticle(partons_[0]->CC(),true)));
    newqbar->set5Momentum(BGFMomenta_[2]);
    ShowerParticlePtr newg   (new_ptr(ShowerParticle(gluon_,false)));
    newg->set5Momentum(BGFMomenta_[0]);
    ShowerParticlePtr newspace(new_ptr(ShowerParticle(partons_[0],false)));
    newspace->set5Momentum(BGFMomenta_[0]-BGFMomenta_[2]);
    HardBranchingPtr spaceBranch(new_ptr(HardBranching(newg,sudakov,HardBranchingPtr(),
						       HardBranching::Incoming)));
    HardBranchingPtr offBranch(new_ptr(HardBranching(newspace,SudakovPtr(),spaceBranch,
						     HardBranching::Incoming)));
    HardBranchingPtr qbar(new_ptr(HardBranching(newqbar,SudakovPtr(),spaceBranch,
						HardBranching::Outgoing)));
    spaceBranch->addChild(offBranch);
    spaceBranch->addChild(qbar);
    HardBranchingPtr outBranch(new_ptr(HardBranching(newq,SudakovPtr(),
						     HardBranchingPtr(),
						     HardBranching::Outgoing)));
    spaceBranchings.push_back(spaceBranch);
    allBranchings.push_back(offBranch);
    allBranchings.push_back(outBranch);
    ColinePtr newin(new_ptr(ColourLine())),newout(new_ptr(ColourLine()));
    newin ->addColoured(newg    ,partons_[0]->id()<0);
    newin ->addColoured(newspace,partons_[0]->id()<0);
    newin ->addColoured(newq    ,partons_[0]->id()<0);
    newout->addColoured(newg    ,partons_[0]->id()>0);
    newout->addColoured(newqbar ,partons_[0]->id()>0);
  }
  HardTreePtr newTree(new_ptr(HardTree(allBranchings,spaceBranchings,
				       ShowerInteraction::QCD)));
  // Set the maximum pt for all other emissions and connect hard and shower tree
  Energy pT = isCompton ? pTCompton_ : pTBGF_;
  // incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    // set maximum pT
    if(QuarkMatcher::Check(cit->first->progenitor()->data()))
      cit->first->maximumpT(pT);
    for(set<HardBranchingPtr>::iterator cjt=newTree->branchings().begin();
	cjt!=newTree->branchings().end();++cjt) {
      if(!(*cjt)->branchingParticle()->isFinalState()&&
	 (*cjt)->branchingParticle()->id()==cit->first->progenitor()->id()) {
	newTree->connect(cit->first->progenitor(),*cjt);
	tPPtr beam =cit->first->original();
	if(!beam->parents().empty()) beam=beam->parents()[0];
	(*cjt)->beam(beam);
	HardBranchingPtr parent=(*cjt)->parent();
	while(parent) {
	  parent->beam(beam);
	  parent=parent->parent();
	};
      }
    }
  }
  // outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    // set maximum pT
    if(QuarkMatcher::Check(cit->first->progenitor()->data()))
      cit->first->maximumpT(pT);
    for(set<HardBranchingPtr>::iterator cjt=newTree->branchings().begin();
	cjt!=newTree->branchings().end();++cjt) {
      if((*cjt)->branchingParticle()->isFinalState()&&
	 (*cjt)->branchingParticle()->id()==cit->first->progenitor()->id()) {
	newTree->connect(cit->first->progenitor(),*cjt);
      }
    }
  }
  // set the evolution partners and scales
  ShowerParticleVector particles;
  for(set<HardBranchingPtr>::iterator cit=newTree->branchings().begin();
      cit!=newTree->branchings().end();++cit) {
    particles.push_back((*cit)->branchingParticle());
  }
  evolver()->showerModel()->partnerFinder()->
    setInitialEvolutionScales(particles,true,ShowerInteraction::QCD,true);
  // Calculate the shower variables:
  evolver()->showerModel()->kinematicsReconstructor()->
    deconstructHardJets(newTree,evolver(),ShowerInteraction::QCD);
  return newTree;
}

void DISHardGenerator::generateCompton() {
  // maximum value of the xT
  double xT = (1.-xB_)/xB_;
  double xTMin = 2.*pTmin_/sqrt(q2_);
  double zp;
  // prefactor
  double a = alphaS_->overestimateValue()*comptonWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.),phi(0.);
  do {
    // intergration variables dxT/xT^3
    xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
    // dzp
    zp = UseRandom::rnd();
    xp = 1./(1.+0.25*sqr(xT)/zp/(1.-zp));
    // check allowed
    if(xp<xB_||xp>1.) continue;
    // other phase-space variables
    phi=UseRandom::rnd()*Constants::twopi;
    // phase-space piece of the weight
    wgt = 8.*(1.-xp)*zp/comptonWeight_;
    // PDF piece of the weight
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1.);
    wgt *= pdf_->xfx(beam_,partons_[0],scale,xB_/xp)/
           pdf_->xfx(beam_,partons_[0],q2_  ,xB_);
    // me piece of the weight
    wgt *= comptonME(xT,xp,zp,phi);
    if(wgt>1.||wgt<0.) generator()->log() << "Compton weight problem " << wgt << "\n";
  }
  while(xT>xTMin&&UseRandom::rnd()>wgt);
  if(xT<=xTMin) {
    pTCompton_=-GeV;
    return;
  }
  // momenta for the configuration
  Energy Q(sqrt(q2_));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		       -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  pTCompton_=0.5*Q*xT;
  ComptonMomenta_.resize(3);
  ComptonMomenta_[0]=p0;
  ComptonMomenta_[1]=p1;
  ComptonMomenta_[2]=p2;
  ComptonISFS_ = zp>xp;
}

void DISHardGenerator::generateBGF() {
  // maximum value of the xT
  double xT = (1.-xB_)/xB_;
  double xTMin = 2.*max(pTmin_,pTCompton_)/sqrt(q2_);
  double zp;
  // prefactor
  double a = alphaS_->overestimateValue()*BGFWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.),phi(0.);
  do {
    // intergration variables dxT/xT^3
    xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
    // dzp
    zp = UseRandom::rnd();
    xp = 1./(1.+0.25*sqr(xT)/zp/(1.-zp));
    // check allowed
    if(xp<xB_||xp>1.) continue;
    // other phase-space variables
    phi=UseRandom::rnd()*Constants::twopi;
    // phase-space piece of the weight
    wgt = 8.*sqr(1.-xp)*zp/BGFWeight_;
    // PDF piece of the weight
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1.);
    wgt *= pdf_->xfx(beam_,gluon_,scale,xB_/xp)/
           pdf_->xfx(beam_,partons_[0],q2_  ,xB_);
    // me piece of the weight
    wgt *= BGFME(xT,xp,zp,phi);
    if(wgt>1.||wgt<0.) generator()->log() << "BGF weight problem " << wgt << "\n";
  }
  while(xT>xTMin&&UseRandom::rnd()>wgt);
  if(xT<=xTMin) {
    pTBGF_=-GeV;
    return;
  }
  // momenta for the configuration
  Energy Q(sqrt(q2_));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		       -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  pTBGF_=0.5*Q*xT;
  BGFMomenta_.resize(3);
  BGFMomenta_[0]=p0;
  BGFMomenta_[1]=p1;
  BGFMomenta_[2]=p2;
}

double DISHardGenerator::comptonME(double xT, double xp, double zp,
				   double phi) {
  // kinematical variables
//   double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
//   double x3 = 2.+x1-x2;
  // matrix element
  vector<double> output(3,0.);
  double cos2 = x2/sqrt(sqr(x2)+sqr(xT));
  double sin2 = xT/sqrt(sqr(x2)+sqr(xT));
  double root = sqrt(sqr(l_)-1.);
  double cphi = cos(phi);
  // compute the r2 term
  double R2= (sqr(cos2)+acoeff_*cos2*(l_-root*sin2*cphi)+sqr(l_-root*sin2*cphi))/
    (1+acoeff_*l_+sqr(l_));
  return 4./3.*alphaS_->ratio(0.25*q2_*sqr(xT))*(1.+R2*sqr(xp)*(sqr(x2)+sqr(xT)));
}

double DISHardGenerator::BGFME(double xT, double xp, double zp,
			       double phi) {
  // kinematical variables
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  // matrix element
  vector<double> output(3,0.);
  double cos2 = x2/sqrt(sqr(x2)+sqr(xT));
  double sin2 = xT/sqrt(sqr(x2)+sqr(xT));
  double cos3 = x3/sqrt(sqr(x3)+sqr(xT));
  double sin3 = xT/sqrt(sqr(x3)+sqr(xT));
  double root = sqrt(sqr(l_)-1.);
  double cphi = cos(phi);
  // compute the r2 term
  double R2= (sqr(cos2)+acoeff_*cos2*(l_-root*sin2*cphi)+sqr(l_-root*sin2*cphi))/
    (1+acoeff_*l_+sqr(l_));
  double R3= (sqr(cos3)-acoeff_*cos3*(l_+root*sin3*cphi)+sqr(l_+root*sin3*cphi))/
    (1+acoeff_*l_+sqr(l_));
  return 0.5*alphaS_->ratio(0.25*q2_*sqr(xT))*sqr(xp)*
    (R2*(sqr(x2)+sqr(xT))+R3*(sqr(x3)+sqr(xT)));
}

double DISHardGenerator::A(tcPDPtr lin, tcPDPtr lout,
			   tcPDPtr qin, tcPDPtr) {
  double output;
  // charged current
  if(lin->id()!=lout->id()) {
    output = 2;
  }
  // neutral current
  else {
    double fact = 0.25*q2_/(q2_+mz2_)/sinW_/cosW_;
    double cvl,cal,cvq,caq;
    if(abs(lin->id())%2==0) {
      cvl = generator()->standardModel()->vnu()*fact+generator()->standardModel()->enu();
      cal = generator()->standardModel()->anu()*fact;
    }
    else {
      cvl = generator()->standardModel()->ve()*fact+generator()->standardModel()->ee();
      cal = generator()->standardModel()->ae()*fact;
    }
    if(abs(qin->id())%2==0) {
      cvq = generator()->standardModel()->vu()*fact+generator()->standardModel()->eu();
      caq = generator()->standardModel()->au()*fact;
    }
    else {
      cvq = generator()->standardModel()->vd()*fact+generator()->standardModel()->ed();
      caq = generator()->standardModel()->ad()*fact;
    }
    output = 8.*cvl*cal*cvq*caq/(sqr(cvl)+sqr(cal))/(sqr(cvq)+sqr(caq));
  }
  if(qin->id()<0) output *= -1.;
  if(lin->id()<0) output *= -1;
  return output;
}
