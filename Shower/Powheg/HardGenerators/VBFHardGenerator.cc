// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VBFHardGenerator class.
//

#include "VBFHardGenerator.h"
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

VBFHardGenerator::VBFHardGenerator() : comptonWeight_(50.), BGFWeight_(150.), 
				       pTmin_(1.*GeV) 
{}

IBPtr VBFHardGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr VBFHardGenerator::fullclone() const {
  return new_ptr(*this);
}

void VBFHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << alphaS_ << ounit(pTmin_,GeV) << comptonWeight_ << BGFWeight_ << gluon_;
}

void VBFHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> alphaS_ >> iunit(pTmin_,GeV) >> comptonWeight_ >> BGFWeight_ >> gluon_;
}

void VBFHardGenerator::doinit() {
  HardestEmissionGenerator::doinit();
  gluon_ = getParticleData(ParticleID::g);
}

ClassDescription<VBFHardGenerator> VBFHardGenerator::initVBFHardGenerator;
// Definition of the static class description member.

void VBFHardGenerator::Init() {

  static ClassDocumentation<VBFHardGenerator> documentation
    ("There is no documentation for the VBFHardGenerator class");

  static Reference<VBFHardGenerator,ShowerAlpha> interfaceShowerAlphaQCD
    ("ShowerAlphaQCD",
     "The object calculating the strong coupling constant",
     &VBFHardGenerator::alphaS_, false, false, true, false, false);

  static Parameter<VBFHardGenerator,Energy> interfacepTMin
    ("pTMin",
     "The minimum pT",
     &VBFHardGenerator::pTmin_, GeV, 1.*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<VBFHardGenerator,double> interfaceComptonWeight
    ("ComptonWeight",
     "Weight for the overestimate ofthe compton channel",
     &VBFHardGenerator::comptonWeight_, 50.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<VBFHardGenerator,double> interfaceBGFWeight
    ("BGFWeight",
     "Weight for the overestimate of the BGF channel",
     &VBFHardGenerator::BGFWeight_, 100.0, 0.0, 1000.0,
     false, false, Interface::limited);

}

bool VBFHardGenerator::canHandle(ShowerTreePtr tree) {
  // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // two outgoing particles
  if(tree->outgoingLines().size()!=3) return false;
  pair<tPPtr,tPPtr> first,second;
  // get the incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(!first.first)  first.first = cit->first->progenitor();
    else             second.first = cit->first->progenitor();
  }
  // incoming quarks or antiquarks
  if(!first.first||!second.first) return false;
  if(abs(first.first->id())>5||abs(second.first->id())>5) return false;
  bool higgs=false;
  // and the outgoing
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
 	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    if(cjt->first->progenitor()->id()==ParticleID::h0) {
      higgs = true;
    }
    else {
      if(abs(cjt->first->progenitor()->id())>5) continue;
      if(cjt->first->progenitor()->colourLine()&&
	 cjt->first->progenitor()->colourLine()==first.first->colourLine()) {
	first.second = cjt->first->progenitor();
	continue;
      }
      if(cjt->first->progenitor()->antiColourLine()&&
	 cjt->first->progenitor()->antiColourLine()==first.first->antiColourLine()) {
	first.second = cjt->first->progenitor();
	continue;
      }
      if(cjt->first->progenitor()->colourLine()&&
	 cjt->first->progenitor()->colourLine()==second.first->colourLine()) {
	second.second = cjt->first->progenitor();
	continue;
      }
      if(cjt->first->progenitor()->antiColourLine()&&
	 cjt->first->progenitor()->antiColourLine()==second.first->antiColourLine()) {
	second.second = cjt->first->progenitor();
	continue;
      }
    }
  }
  if(!first.second||!second.second||!higgs) return false;
  // can handle it
  return true;
}

HardTreePtr VBFHardGenerator::generateHardest(ShowerTreePtr tree) {
  pair<    tShowerParticlePtr,    tShowerParticlePtr> first,second;
  pair<tcBeamPtr,tcBeamPtr> beams;
  pair<tPPtr,tPPtr> hadrons;
  tShowerParticlePtr higgs;
  // get the incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(!first.first) {
      first.first    = cit->first->progenitor();
      beams.first    = cit->first->beam();
      hadrons.first  = cit->first->original()->parents()[0];
    }
    else {
      second.first   = cit->first->progenitor();
      beams.second   = cit->first->beam();
      hadrons.second = cit->first->original()->parents()[0];
    }
  }
  // and the outgoing
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
 	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    if(cjt->first->progenitor()->id()==ParticleID::h0) {
      higgs = cjt->first->progenitor();
    }
    else {
      if(abs(cjt->first->progenitor()->id())>5) continue;
      if(cjt->first->progenitor()->colourLine()&&
	 cjt->first->progenitor()->colourLine()==first.first->colourLine()) {
	first.second = cjt->first->progenitor();
	continue;
      }
      if(cjt->first->progenitor()->antiColourLine()&&
	 cjt->first->progenitor()->antiColourLine()==first.first->antiColourLine()) {
	first.second = cjt->first->progenitor();
	continue;
      }
      if(cjt->first->progenitor()->colourLine()&&
	 cjt->first->progenitor()->colourLine()==second.first->colourLine()) {
	second.second = cjt->first->progenitor();
	continue;
      }
      if(cjt->first->progenitor()->antiColourLine()&&
	 cjt->first->progenitor()->antiColourLine()==second.first->antiColourLine()) {
	second.second = cjt->first->progenitor();
	continue;
      }
    }
  }
  // first system will emit random swap so both with equal probability
  if(UseRandom::rnd()<0.5) {
    swap(first,second);
    swap(beams.first,beams.second);
  }
  // check beam, all particles
  assert(beams.first  && higgs &&
	 first .first &&  first.second && 
	 second.first && second.second);
  // beam and pdf
  beam_ = beams.first;
  pdf_  = beam_->pdf();
  assert(beam_&&pdf_);
  // Particle data objects
  partons_[0] =  first. first->dataPtr();
  partons_[1] =  first.second->dataPtr();
  partons_[2] = second. first->dataPtr();
  partons_[3] = second.second->dataPtr();
  higgs->dataPtr();
  // extract the born variables
  q_ = first.second->momentum()-first.first->momentum();
  q2_ = -q_.m2();
  xB_ = first.first->x();
  Lorentz5Momentum pb     = first.first->momentum();
  Lorentz5Momentum pbasis = hadrons.first->momentum();
  pbasis.setMass(0.*GeV);
  pbasis.rescaleRho();
  Axis axis(q_.vect().unit());
  double sinth(sqr(axis.x())+sqr(axis.y()));
  rot_ = LorentzRotation();
  if(axis.perp2()>1e-20) {
    rot_.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    rot_.rotateX(Constants::pi);
  }
  if(abs(1.-q_.e()/q_.vect().mag())>1e-6) rot_.boostZ( q_.e()/q_.vect().mag());
  pb *= rot_;
  if(pb.perp2()/GeV2>1e-20) {
    Boost trans = -1./pb.e()*pb.vect();
    trans.setZ(0.);
    rot_.boost(trans);
  }
  // momenta of the particles
  phiggs_     = rot_*higgs->momentum();
  pother_ [0] = rot_*second. first->momentum();
  pother_ [1] = rot_*second.second->momentum();
  psystem_[0] = rot_* first. first->momentum();
  psystem_[1] = rot_* first.second->momentum();
  q_ *= rot_;
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
				 << "VBFHardGenerator::generateHardest()" 
				 << Exception::runerror;
  // add the non emitting particles
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  spaceBranchings.push_back(new_ptr(HardBranching(second.first,SudakovPtr(),
						  HardBranchingPtr(),
						  HardBranching::Incoming)));
  allBranchings.push_back(spaceBranchings.back());
  allBranchings.push_back(new_ptr(HardBranching(second.second,SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  allBranchings.push_back(new_ptr(HardBranching(higgs,SudakovPtr(),
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
      ColinePtr newline(new_ptr(ColourLine()));
      newline->addColoured(newspace,newspace->dataPtr()->iColour()!=PDT::Colour3);
      newline->addColoured(newqout ,newspace->dataPtr()->iColour()!=PDT::Colour3);
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

      ColinePtr newline(new_ptr(ColourLine()));
      newline->addColoured(newqin ,newqin->dataPtr()->iColour()!=PDT::Colour3);
      newline->addColoured(newtime,newqin->dataPtr()->iColour()!=PDT::Colour3);
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

    ColinePtr newline(new_ptr(ColourLine()));
    newline->addColoured(newspace,newspace->dataPtr()->iColour()!=PDT::Colour3);
    newline->addColoured(newq    ,newspace->dataPtr()->iColour()!=PDT::Colour3);
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
    set<HardBranchingPtr>::iterator cjt=newTree->branchings().begin();
    if(cit->first->progenitor()==first.first) {
      ++cjt;++cjt;++cjt;
    }
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
  for(set<HardBranchingPtr>::iterator cjt=newTree->branchings().begin();
      cjt!=newTree->branchings().end();++cjt) {
    if(cjt==newTree->branchings().begin()) {
      (**cjt).showerMomentum((**cjt).branchingParticle()->momentum());
      ++cjt;
      (**cjt).showerMomentum((**cjt).branchingParticle()->momentum());
      ++cjt;
      (**cjt).showerMomentum((**cjt).branchingParticle()->momentum());
      ++cjt;
    }
    // incoming
    if((**cjt).status()==HardBranching::Incoming) {
      first. first->set5Momentum((**cjt).showerMomentum());
    }
    // outgoing
    else {
      first.second->set5Momentum((**cjt).showerMomentum());
    }
  }
  return newTree;
}

void VBFHardGenerator::generateCompton() {
  // maximum value of the xT
  double xT = sqrt((1.-xB_)/xB_);
  double xTMin = 2.*pTmin_/sqrt(q2_);
  double zp;
  // prefactor
  double a = alphaS_->overestimateValue()*comptonWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.),phi(0.);
  do {
    wgt = 0.;
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

void VBFHardGenerator::generateBGF() {
  // maximum value of the xT
  double xT = (1.-xB_)/xB_;
  double xTMin = 2.*max(pTmin_,pTCompton_)/sqrt(q2_);
  double zp;
  // prefactor
  double a = alphaS_->overestimateValue()*BGFWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.),phi(0.);
  do {
    wgt = 0.;
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

double VBFHardGenerator::comptonME(double xT, double xp, double zp,
				   double phi) {
  // scale and prefactors
  Energy2 mu2 = q2_;
  double aS = 2.*alphaS_->value(mu2);
  double CFfact = 4./3.*aS/Constants::twopi;
  double sin2W = generator()->standardModel()->sin2ThetaW();
  Energy Q(sqrt(q2_));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  //set NLO momenta
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		       -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  Lorentz5Momentum qnlo = p2+p1-p0; 
  Energy2 q2nlo = qnlo.m2();
  // Breit frame variables
  Lorentz5Momentum r1 = -p0/x1;
  Lorentz5Momentum r2 =  p1/x2;
  // electroweak parameters
  Energy2 mz2 = sqr(getParticleData(ParticleID::Z0)->mass());
  Energy Gamma(getParticleData(ParticleID::Z0)->width());
  Energy4 MG2(mz2*sqr(Gamma));
  double c0L,c1L,c0R,c1R;
  if(abs(partons_[0]->id())%2==0) {
    c0L = 0.5*generator()->standardModel()->au()-
	  sin2W*generator()->standardModel()->eu();
    c0R =-sin2W*generator()->standardModel()->eu();
  }
  else {
    c0L = 0.5*generator()->standardModel()->ad()-
	  sin2W*generator()->standardModel()->ed();
    c0R =-sin2W*generator()->standardModel()->ed();
  }

  if(abs(partons_[2]->id())%2==0) {
    c1L = 0.5*generator()->standardModel()->au()-
	  sin2W*generator()->standardModel()->eu();
    c1R =-sin2W*generator()->standardModel()->eu();
  }
  else {
    c1L = 0.5*generator()->standardModel()->ad()-
	  sin2W*generator()->standardModel()->ed();
    c1R =-sin2W*generator()->standardModel()->ed();
  }
  // Matrix element variables
  double G1 = sqr(c0L*c1L)+sqr(c0R*c1R);
  double G2 = sqr(c0L*c1R)+sqr(c0R*c1L);
  InvEnergy4 Dlo = 1./(sqr(q2_  -mz2)+MG2); 
  InvEnergy4 Dbf = 1./(sqr(q2nlo-mz2)+MG2);

  Energy4 term1 = G1*(r1*pother_[0])*((qnlo+r1)*pother_[1])+
                  G2*(r1*pother_[1])*((qnlo+r1)*pother_[0]);

  Energy4 term2 = G1*((r2-qnlo)*pother_[0])*(r2*pother_[1])+
                  G2*((r2-qnlo)*pother_[1])*(r2*pother_[0]);

  Energy4 commfact2 = G1*(psystem_[0]*pother_[0])*
                         (psystem_[1]*pother_[1])+
                      G2*(psystem_[0]*pother_[1])*
                         (psystem_[1]*pother_[0]);

  return 1./3.*CFfact*1./((1.-xp)*(1.-zp))*sqr(Dbf/Dlo)
                    /commfact2*(term1+sqr(xp-1.+zp)*term2); 
}

double VBFHardGenerator::BGFME(double xT, double xp, double zp,
			       double phi) {
  // scale and prefactors
  Energy2 mu2 = q2_;
  double aS = 2.*alphaS_->value(mu2);
  double TRfact = 1./2.*aS/Constants::twopi;
  double sin2W = generator()->standardModel()->sin2ThetaW();
  Energy Q(sqrt(q2_));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  // Set NLO momenta
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		       -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  Lorentz5Momentum qnlo = p2+p1-p0; 
  Energy2 q2nlo = qnlo.m2();
  // Breit frame variables
  Lorentz5Momentum r1 = -p0/x1;
  Lorentz5Momentum r2 =  p1/x2;
  // electroweak parameters
  Energy2 mz2 = sqr(getParticleData(ParticleID::Z0)->mass());
  Energy Gamma(getParticleData(ParticleID::Z0)->width());
  Energy4 MG2(mz2*sqr(Gamma));
  double c0L,c1L,c0R,c1R;
  if(abs(partons_[0]->id())%2==0) {
    c0L = 0.5*generator()->standardModel()->au()-
	  sin2W*generator()->standardModel()->eu();
    c0R =-sin2W*generator()->standardModel()->eu();
  }
  else {
    c0L = 0.5*generator()->standardModel()->ad()-
	  sin2W*generator()->standardModel()->ed();
    c0R =-sin2W*generator()->standardModel()->ed();
  }

  if(abs(partons_[2]->id())%2==0) {
    c1L = 0.5*generator()->standardModel()->au()-
	  sin2W*generator()->standardModel()->eu();
    c1R =-sin2W*generator()->standardModel()->eu();
  }
  else {
    c1L = 0.5*generator()->standardModel()->ad()-
	  sin2W*generator()->standardModel()->ed();
    c1R =-sin2W*generator()->standardModel()->ed();
  }
  // Matrix element variables
  double G1 = sqr(c0L*c1L)+sqr(c0R*c1R);
  double G2 = sqr(c0L*c1R)+sqr(c0R*c1L);
  InvEnergy4 Dlo = 1./(sqr(q2_  -mz2)+MG2); 
  InvEnergy4 Dbf = 1./(sqr(q2nlo-mz2)+MG2);

  Energy4 term1 = G1*(r1*pother_[0])*((qnlo+r1)*pother_[1])+
                  G2*(r1*pother_[1])*((qnlo+r1)*pother_[0]);

  Energy4 term2 = G1*((r2-qnlo)*pother_[0])*(r2*pother_[1])+
                  G2*((r2-qnlo)*pother_[1])*(r2*pother_[0]);

  Energy4 commfact2 = G1*(psystem_[0]*pother_[0])*
                         (psystem_[1]*pother_[1])+
                      G2*(psystem_[0]*pother_[1])*
                         (psystem_[1]*pother_[0]);

  return  1./3.*TRfact*1./((1.-zp)*sqr(xp)*
          (2+xp-zp))* sqr(Dbf/Dlo)/commfact2*
          (term1+sqr(xp-1.+zp)*term2);
}

