// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGammaHardGenerator class.
//

#include "GammaGammaHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Cuts/Cuts.h"
#include <numeric>

using namespace Herwig;

GammaGammaHardGenerator::GammaGammaHardGenerator()  
  : pTmin_(2.*GeV), qqgFactor_(10.), qgqFactor_(12.), qbargqbarFactor_(12.),
    power_(1.), minPhotonpT_(-GeV)
{}

IBPtr GammaGammaHardGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGammaHardGenerator::fullclone() const {
  return new_ptr(*this);
}

void GammaGammaHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << alphaS_ << gluon_ << photon_ << ounit(pTmin_,GeV) << power_
     << FFPvertex_ << FFGvertex_ << qqgFactor_ << qgqFactor_ << qbargqbarFactor_;
}

void GammaGammaHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> alphaS_ >> gluon_ >> photon_ >> iunit(pTmin_,GeV) >> power_
     >> FFPvertex_ >> FFGvertex_ >> qqgFactor_ >> qgqFactor_ >> qbargqbarFactor_;
}

// Definition of the static class description member.
ClassDescription<GammaGammaHardGenerator> 
GammaGammaHardGenerator::initGammaGammaHardGenerator;

void GammaGammaHardGenerator::Init() {

  static ClassDocumentation<GammaGammaHardGenerator> documentation
    ("The GammaGammaHardGenerator class implements the generation of hard"
     " QCD radiation in q qbar -> gamma gamma events");

  static Reference<GammaGammaHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &GammaGammaHardGenerator::alphaS_, false, false, true, false, false);

  static Parameter<GammaGammaHardGenerator, Energy> interfacePtMin
    ("minPt",
     "The pt cut on hardest emision generation",
     &GammaGammaHardGenerator::pTmin_, GeV, 2.*GeV, ZERO, 100000.0*GeV,
     false, false, Interface::limited);

  static Parameter<GammaGammaHardGenerator,double> interfaceSamplingPower
    ("SamplingPower",
     "Power of pT for the overestimate in the Sudakov",
     &GammaGammaHardGenerator::power_, 1., 0.0, 2.0,
     false, false, Interface::limited);

  static Parameter<GammaGammaHardGenerator,double> interfaceQQBarFactor
    ("QQBarFactor",
     "Factor for the overestimate of the q qbar channel",
     &GammaGammaHardGenerator::qqgFactor_, 1.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<GammaGammaHardGenerator,double> interfaceQGFactor
    ("QGFactor",
     "Factor for the overestimate of the q g channel",
     &GammaGammaHardGenerator::qgqFactor_, 1.0,  0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<GammaGammaHardGenerator,double> interfacegQBarFactor
    ("gQBarFactor",
     "Factor for the overestimate of the g qbar channel",
     &GammaGammaHardGenerator::qbargqbarFactor_, 1.0, 0.0, 1000.0,
     false, false, Interface::limited);

}

void GammaGammaHardGenerator::doinit() {
  HardestEmissionGenerator::doinit();
  gluon_  = getParticleData(ParticleID::g);
  photon_ = getParticleData(ParticleID::gamma);
  // get a pointer to the standard model object in the run
  static const tcHwSMPtr hwsm
    = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if (!hwsm) throw InitException() << "hwsm pointer is null in"
				   << " MEPP2GammaGamma::doinit()"
				   << Exception::abortnow;
  // get pointers to all required Vertex objects
  FFPvertex_ = hwsm->vertexFFP();
  FFGvertex_ = hwsm->vertexFFG();
}


void GammaGammaHardGenerator::doinitrun() {
  HardestEmissionGenerator::doinitrun();
  minPhotonpT_ = generator()->eventHandler()->cuts()->minKT(photon_);
  minPhotonEta_ = generator()->eventHandler()->cuts()->minEta(photon_);
  maxPhotonEta_ = generator()->eventHandler()->cuts()->maxEta(photon_);
}

bool GammaGammaHardGenerator::canHandle(ShowerTreePtr tree) {
  // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // should be a quark and an antiquark
  unsigned int ix(0);
  ShowerParticlePtr part[2];
  // extract the incoming particles
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    part[ix]=cit->first->progenitor();
    ++ix;
  }
  // check incoming quark and antiquark
  if(!(QuarkMatcher::Check(part[0]->data())&&QuarkMatcher::Check(part[1]->data())&&
       double(part[0]->id())/double(part[1]->id())<0.)) return false;
  // two outgoing particles
  if(tree->outgoingLines().size()!=2) return false;
  // extract the outgoing particles
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    if(cjt->first->progenitor()->id()!=ParticleID::gamma) return false;
  }
  return true;
}

HardTreePtr GammaGammaHardGenerator::generateHardest(ShowerTreePtr tree) {
  // get the particles to be showered
  beams_.clear();
  partons_.clear();
  // find the incoming particles
  ShowerParticleVector incoming;
  vector<ShowerProgenitorPtr> particlesToShower;
  //progenitor particles are produced in z direction.
  for( map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	 cit = tree->incomingLines().begin(); 
       cit != tree->incomingLines().end(); ++cit ) {
    incoming.push_back( cit->first->progenitor() );
    beams_.push_back( cit->first->beam() );
    partons_.push_back( cit->first->progenitor()->dataPtr() );
    particlesToShower.push_back( cit->first );
  }
  // find the outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
 	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    particlesToShower.push_back(cjt->first);
  }
  for(unsigned int ix=0;ix<particlesToShower.size();++ix)
  // CMS energy of the hadron collision
  s_  = generator()->currentEvent()->primaryCollision()->m2();
  rs_ = sqrt(s_);
  // Borm momentum fractions
  Lorentz5Momentum psystem(particlesToShower[2]->progenitor()->momentum()+
			   particlesToShower[3]->progenitor()->momentum());
  systemMass_ = psystem.m();
  double systemRapidity = psystem.rapidity();
  x_[0] = systemMass_*exp( systemRapidity)/sqrt(s_);          
  x_[1] = systemMass_*exp(-systemRapidity)/sqrt(s_);
  // leading order momenta
  for(unsigned int ix=0;ix<4;++ix)
    loMomenta_[ix] = particlesToShower[ix]->progenitor()->momentum();
  // generate the new configuration
  pTqqbar_[0] = -GeV;
  pTqqbar_[1] = -GeV;
  pTqg_       = -GeV;
  pTgqbar_    = -GeV;
  // generate a q qbar -> gamma gamma g    configuration
  generateQQbarG();
  // generate a q g    -> gamma gamma q    configuration
  generateQGQ();
  // generate a g qbar -> gamma gamma qbar configuration
  generateQbarGQbar();
  // return if no emission
  if(pTqqbar_[0]<ZERO&&pTqqbar_[1]<ZERO&&pTqg_<ZERO&&pTgqbar_<ZERO) {
    for(unsigned int ix=0;ix<particlesToShower.size();++ix)
      particlesToShower[ix]->maximumpT(pTmin_);
    return HardTreePtr();
  }
  // select the type of emission
  int emissionType=0, iemit(0);
  if      (pTqqbar_[0]>pTqqbar_[1] && pTqqbar_[0]>pTqg_    && 
	   pTqqbar_[0]>pTgqbar_ ) {
    emissionType = 1;
    iemit = 0;
  }
  else if (pTqqbar_[1]>pTqg_       && pTqqbar_[1]>pTgqbar_ ) {
    emissionType = 1;
    iemit = 1;
  }
  else if (pTqg_>pTgqbar_    ) {
    emissionType = 2;
  }
  else {
    emissionType = 3;
  }
  ShowerParticleVector newparticles;
  Lorentz5Momentum pphoton[2];
  Energy pTEmit;
  // q qbar -> V gamma g
  if(emissionType==1) {
    newparticles.push_back(new_ptr(ShowerParticle(partons_[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(partons_[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gluon_           , true)));
    if(partons_[0]->id()>0) {
      newparticles[0]->set5Momentum(pQqqbar_   [iemit]);
      newparticles[1]->set5Momentum(pQbarqqbar_[iemit]);
    }
    else {
      newparticles[0]->set5Momentum(pQbarqqbar_[iemit]);
      newparticles[1]->set5Momentum(pQqqbar_   [iemit]);
    }
    newparticles[2]->set5Momentum(pGqqbar_[iemit]);
    pphoton[0] = pGammaqqbar_[iemit][0];
    pphoton[1] = pGammaqqbar_[iemit][1];
    pTEmit = pTqqbar_[iemit];
  }
  // q g    -> q V gamma
  else if(emissionType==2) {
    if(partons_[0]->id()>0) {
      newparticles.push_back(new_ptr(ShowerParticle(partons_[0]      ,false)));
      newparticles.push_back(new_ptr(ShowerParticle(gluon_           ,false)));
      newparticles.push_back(new_ptr(ShowerParticle(partons_[1]->CC(), true)));
      newparticles[0]->set5Momentum(pQinqg_);
      newparticles[1]->set5Momentum(pGqg_);
      iemit=1;
    }
    else {
      newparticles.push_back(new_ptr(ShowerParticle(gluon_           ,false)));
      newparticles.push_back(new_ptr(ShowerParticle(partons_[1]      ,false)));
      newparticles.push_back(new_ptr(ShowerParticle(partons_[0]->CC(), true)));
      newparticles[1]->set5Momentum(pQinqg_);
      newparticles[0]->set5Momentum(pGqg_);
      iemit=0;
    }
    newparticles[2]->set5Momentum(pQoutqg_);
    pphoton[0] = pGammaqg_[0];
    pphoton[1] = pGammaqg_[1];
    pTEmit = pTqg_;
  }
  // g qbar -> qbar V gamma
  else {
    iemit=0;
    if(partons_[0]->id()>0) {
      generator()->log() << "testing did A\n";
      newparticles.push_back(new_ptr(ShowerParticle(gluon_           ,false)));
      newparticles.push_back(new_ptr(ShowerParticle(partons_[1]      ,false)));
      newparticles.push_back(new_ptr(ShowerParticle(partons_[0]->CC(), true)));
      newparticles[0]->set5Momentum(pGgqbar_);
      newparticles[1]->set5Momentum(pQingqbar_);
      iemit=0;
    }
    else {
      generator()->log() << "testing did B\n";
      newparticles.push_back(new_ptr(ShowerParticle(partons_[0]      ,false)));
      newparticles.push_back(new_ptr(ShowerParticle(gluon_           ,false)));
      newparticles.push_back(new_ptr(ShowerParticle(partons_[1]->CC(), true)));
      newparticles[1]->set5Momentum(pGgqbar_);
      newparticles[0]->set5Momentum(pQingqbar_);
      iemit=1;
    }
    newparticles[2]->set5Momentum(pQoutgqbar_);
    pphoton[0] = pGammagqbar_[0];
    pphoton[1] = pGammagqbar_[1];
    pTEmit  = pTgqbar_;
  }
  // create the photons
  for(unsigned int ix=0;ix<2;++ix) {
    newparticles.push_back(new_ptr(ShowerParticle(photon_,true)));
    newparticles.back()->set5Momentum(pphoton[ix]);
  }
  Lorentz5Momentum poff=newparticles[iemit]->momentum()-newparticles[2]->momentum();
  poff.rescaleMass();
  newparticles.push_back(new_ptr(ShowerParticle(partons_[iemit],false)));
  newparticles.back()->set5Momentum(poff);
  // find the sudakov for the branching
  BranchingList branchings=evolver()->splittingGenerator()->initialStateBranchings();
  long index = abs(partons_[iemit]->id());
  IdList br(3);
  // types of particle in the branching
  br[0]=newparticles[iemit]->id();
  br[1]=newparticles[  5  ]->id();
  br[2]=newparticles[  2  ]->id();
  if(emissionType==1) {
    br[0]=abs(br[0]);
    br[1]=abs(br[1]);
  }
  else if(emissionType==2) {
    br[1]=-br[1];
    br[2]=-br[2];
  }
  SudakovPtr sudakov;
  for(BranchingList::const_iterator cit = branchings.lower_bound(index); 
      cit != branchings.upper_bound(index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0]==br[0]&&ids[1]==br[1]&&ids[2]==br[2]) {
      sudakov=cit->second.first;
      break;
    }
  }
  if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
				 << "GammaGammaHardGenerator::generateHardest()" 
				 << Exception::runerror;
  vector<HardBranchingPtr> nasonin,nasonhard;
  // create the branchings for the incoming particles
  nasonin.push_back(new_ptr(HardBranching(newparticles[0],
					  iemit==0 ? sudakov : SudakovPtr(),
					  HardBranchingPtr(),HardBranching::Incoming)));
  nasonin.push_back(new_ptr(HardBranching(newparticles[1],
					  iemit==1 ? sudakov : SudakovPtr(),
					  HardBranchingPtr(),HardBranching::Incoming)));
  // create the branching for the emitted jet
  nasonin[iemit]->addChild(new_ptr(HardBranching(newparticles[2],SudakovPtr(),
						 nasonin[iemit],
						 HardBranching::Outgoing)));
  // intermediate IS particle
  nasonhard.push_back(new_ptr(HardBranching(newparticles[5],SudakovPtr(),
					    nasonin[iemit],HardBranching::Incoming)));
  nasonin[iemit]->addChild(nasonhard.back());
  // set the colour partners
  nasonhard.back()->colourPartner(nasonin[iemit==0 ? 1 : 0]);
  nasonin[iemit==0 ? 1 : 0]->colourPartner(nasonhard.back());
  // add other particle
  nasonhard.push_back(nasonin[iemit==0 ? 1 : 0]);
  // outgoing photons
  nasonhard.push_back(new_ptr(HardBranching(newparticles[3],SudakovPtr(),
					    HardBranchingPtr(),HardBranching::Outgoing)));
  nasonhard.push_back(new_ptr(HardBranching(newparticles[4],SudakovPtr(),
					    HardBranchingPtr(),HardBranching::Outgoing)));
  // make the tree
  HardTreePtr hardTree=new_ptr(HardTree(nasonhard,nasonin,ShowerInteraction::QCD));
  // connect the ShowerParticles with the branchings
  // and set the maximum pt for the radiation
  set<HardBranchingPtr> hard=hardTree->branchings();
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if( pTEmit < pTmin_ ) particlesToShower[ix]->maximumpT(pTmin_);
    else particlesToShower[ix]->maximumpT(pTmin_);
    for(set<HardBranchingPtr>::const_iterator mit=hard.begin();
	mit!=hard.end();++mit) {
      if(particlesToShower[ix]->progenitor()->id()==(*mit)->branchingParticle()->id()&&
	 (( (*mit)->status()==HardBranching::Incoming &&
	    !particlesToShower[ix]->progenitor()->isFinalState())||
	  ( (*mit)->status()==HardBranching::Outgoing&&
	    particlesToShower[ix]->progenitor()->isFinalState()))) {
	hardTree->connect(particlesToShower[ix]->progenitor(),*mit);
	if((*mit)->status()==HardBranching::Incoming) {
	  (*mit)->beam(particlesToShower[ix]->original()->parents()[0]);
	}
	HardBranchingPtr parent=(*mit)->parent();
	while(parent) {
	  parent->beam(particlesToShower[ix]->original()->parents()[0]);
	  parent=parent->parent();
	};
	hard.erase(mit);
	break;
      }
    }
  }
  ColinePtr newline=new_ptr(ColourLine());
  for(set<HardBranchingPtr>::const_iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3)
      newline->addColoured((**cit).branchingParticle());
    else if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3bar)
      newline->addAntiColoured((**cit).branchingParticle());
  }
  ShowerParticleVector particles;
  for(set<HardBranchingPtr>::iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    particles.push_back((*cit)->branchingParticle());
  }
  evolver()->showerModel()->partnerFinder()->
     setInitialEvolutionScales(particles,true,ShowerInteraction::QCD,true);
  // calculate the shower variables
  evolver()->showerModel()->kinematicsReconstructor()->
    deconstructHardJets(hardTree,evolver(),ShowerInteraction::QCD);
  vector<bool> matched(particlesToShower.size(),false);
  for(set<HardBranchingPtr>::const_iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      if(matched[ix]) continue;
      if((((**cit).status()==HardBranching::Incoming && 
	   !particlesToShower[ix]->progenitor()->isFinalState())||
	  ((**cit).status()==HardBranching::Outgoing &&
	   particlesToShower[ix]->progenitor()->isFinalState()))&&
	 particlesToShower[ix]->progenitor()->id()==(**cit).branchingParticle()->id()) {
	particlesToShower[ix]->progenitor()->set5Momentum((**cit).showerMomentum());
	matched[ix] = true;
	break;
      }
    }
  }
  return hardTree;
}

void GammaGammaHardGenerator::generateQQbarG() {
  // PDFs for the LO cross section
  double pdf[4];
  pdf[0]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],sqr(systemMass_),x_[0])/x_[0];
  pdf[1]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],sqr(systemMass_),x_[1])/x_[1];
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0;
  // loop over the emitting particles
  for(unsigned int ix=0;ix<2;++ix) {
    // storage of pt and rapidity of the gluon
    Energy pT = 0.5*rs_;
    double yj,phi,x1,x2;
    // prefactor for veto algorithm
    double c = alphaS_->overestimateValue()/Constants::twopi*
      qqgFactor_*(maxyj-minyj);
    if(power_!=1.) c /= power_-1.;
    double wgt(0.);
    do {
      // generate pT
      if(power_!=1.) {
	pT = rs_*pow(pow(double(pT/rs_),1.-power_)
		     -log(UseRandom::rnd())/c,1./(1.-power_));
      }
      else {
	pT *= pow(UseRandom::rnd(),1./c);
      }
      // generate rapidity of the jet
      yj = UseRandom::rnd()*(maxyj-minyj)+ minyj;
      // generate phi
      phi = UseRandom::rnd()*Constants::twopi;
      // the CS variables
      double z,vt;
      if(ix==0) {
	z  = (1.-pT/rs_/x_[1]*exp(-yj))/(1.+pT/rs_*exp( yj)/x_[0]);
	vt = pT/rs_/x_[1]*exp(-yj);
	x1 = x_[0]/z;
	x2 = x_[1];
	if(z<x_[0]||z>1.||vt>1.-z||vt<0) continue;
      }
      else {
	z = (1.-pT/rs_/x_[0]*exp( yj))/(1.+pT/rs_*exp(-yj)/x_[1]);
	vt = pT/rs_/x_[0]*exp( yj);
	x1 = x_[0];
	x2 = x_[1]/z;
	if(z<x_[1]||z>1.||vt>1.-z||vt<0) continue;
      }
      // momentum of the gluon
      pGqqbar_[ix] = Lorentz5Momentum (pT *cos(phi),pT *sin(phi),
				       pT *sinh(yj),pT *cosh(yj),ZERO);
      // momenta of the incoming quarks
      if(partons_[0]->id()>0) {
	pQqqbar_   [ix]  = Lorentz5Momentum (ZERO,ZERO,
					     x1*0.5*rs_,x1*0.5*rs_,ZERO);
	pQbarqqbar_[ix]  = Lorentz5Momentum (ZERO,ZERO,
					     -x2*0.5*rs_,x2*0.5*rs_,ZERO);
      }
      else {
	pQbarqqbar_[ix]  = Lorentz5Momentum (ZERO,ZERO,
					     x1*0.5*rs_,x1*0.5*rs_,ZERO);
	pQqqbar_   [ix]  = Lorentz5Momentum (ZERO,ZERO,
					     -x2*0.5*rs_,x2*0.5*rs_,ZERO);
      }
      // momenta of the photons
      Lorentz5Momentum Kt = loMomenta_[0]+loMomenta_[1];
      Lorentz5Momentum K  = pQqqbar_   [ix]+pQbarqqbar_[ix]-pGqqbar_[ix];
      Lorentz5Momentum Ksum = K+Kt;
      Energy2 K2    = K   .m2();
      Energy2 Ksum2 = Ksum.m2();
      for(unsigned int iz=0;iz<2;++iz) { 
	pGammaqqbar_[ix][iz] = loMomenta_[2+iz]
	  -2.*(loMomenta_[2+iz]*Ksum)/Ksum2*Ksum+2.*(Kt*loMomenta_[2+iz])/K2*K;
      }
      // enforce LO photon cuts to avoid QED singularities
      if(pGammaqqbar_[ix][0].perp() < minPhotonpT_  ||
	 pGammaqqbar_[ix][1].perp() < minPhotonpT_  ||
	 pGammaqqbar_[ix][0].eta()  < minPhotonEta_ ||
	 pGammaqqbar_[ix][1].eta()  < minPhotonEta_ ||
	 pGammaqqbar_[ix][0].eta()  > maxPhotonEta_ ||
	 pGammaqqbar_[ix][1].eta()  > maxPhotonEta_ ) {
	wgt=0.;
	continue;
      }
      // now for the weight
      // phase space/over estimate
      wgt  = 2./(1.-vt)*pow(pT/rs_,power_+1)/qqgFactor_;
      // pdf bit
      Energy2 scale = sqr(systemMass_)+sqr(pT);
      pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x1)/x1;
      pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,x2)/x2;
      if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
	wgt=0.;
	continue;
      }
      wgt *= pdf[2]*pdf[3]/pdf[0]/pdf[1];
      // alpha S
      wgt *= alphaS_->ratio(sqr(pT));
      // final bit
      wgt *= QQbarGRatio(ix);
      if(wgt>1.) {
	cerr << "testing weight QQBar " << wgt << " " << pT/GeV << " " << yj << "\n";
	cerr << "testing photon momenta " << pGammaqqbar_[ix][0]/GeV << "\n";
	cerr << "testing photon momenta " << pGammaqqbar_[ix][1]/GeV << "\n";
      }
    }
    while(UseRandom::rnd()>wgt&&pT>pTmin_);
    if(pT<pTmin_) pT=-GeV;
    pTqqbar_[ix] = pT;
  }
}

double GammaGammaHardGenerator::QQbarGRatio(unsigned int order) {
  using namespace ThePEG::Helicity;
  double sum(0.);
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qbarin;
  vector<VectorWaveFunction> pout1,pout2,gout;
  SpinorWaveFunction q_in;
  SpinorBarWaveFunction qbar_in;
  if(partons_[0]->id()>0) {
    q_in    = SpinorWaveFunction   (pQqqbar_   [order],partons_[0],incoming);
    qbar_in = SpinorBarWaveFunction(pQbarqqbar_[order],partons_[1],incoming);
  }
  else {
    q_in    = SpinorWaveFunction   (pQqqbar_   [order],partons_[1],incoming);
    qbar_in = SpinorBarWaveFunction(pQbarqqbar_[order],partons_[0],incoming);
  }
  VectorWaveFunction     p_out1(pGammaqqbar_[order][0],photon_,outgoing);
  VectorWaveFunction     p_out2(pGammaqqbar_[order][1],photon_,outgoing);
  VectorWaveFunction     g_out (pGqqbar_[order]       ,gluon_ ,outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q_in.reset(ix);               qin.push_back(q_in   );
    qbar_in.reset(ix);         qbarin.push_back(qbar_in);
    g_out.reset(2*ix);           gout.push_back(g_out  );
    p_out1.reset(2*ix);         pout1.push_back(p_out1 );
    p_out2.reset(2*ix);         pout2.push_back(p_out2 );
  }
  vector<Complex> diag(6);
  Energy2 scale = sqr(systemMass_);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int phel1=0;phel1<2;++phel1) {
	for(unsigned int phel2=0;phel2<2;++phel2) {
	  for(unsigned int ghel=0;ghel<2;++ghel) {
 	    // first diagram
	    SpinorWaveFunction inters1 = 
	      FFPvertex_->evaluate(ZERO,5,qin[ihel1].particle()->CC(),
				   qin[ihel1],pout1[phel1]);
	    SpinorBarWaveFunction inters2 = 
	      FFPvertex_->evaluate(ZERO,5,qbarin[ihel2].particle()->CC(),
				   qbarin[ihel2],pout2[phel2]);
	    diag[0] = FFGvertex_->evaluate(scale,inters1,inters2,gout[ghel]);
	    // second diagram
	    SpinorWaveFunction inters3 = 
	      FFGvertex_->evaluate(scale,5,qin[ihel1].particle()->CC(),
				   qin[ihel1],gout[ghel]);
	    SpinorBarWaveFunction inters4 = 
	      FFPvertex_->evaluate(ZERO,5,qbarin[ihel2].particle()->CC(),
				   qbarin[ihel2],pout1[phel1]);
	    diag[1] = FFPvertex_->evaluate(ZERO,inters3,inters4,pout2[phel2]);
	    // fourth diagram
	    diag[2] = FFPvertex_->evaluate(ZERO,inters3,inters2,pout1[phel1]);
	    // fifth diagram
	    SpinorBarWaveFunction inters5 = 
	      FFGvertex_->evaluate(scale,5,qbarin[ihel2].particle()->CC(),
				   qbarin[ihel2],gout[ghel]);
	    diag[3] = 
	      FFPvertex_->evaluate(ZERO,inters1,inters5,pout2[phel2]);
	    // sixth diagram
	    SpinorWaveFunction inters6 = 
	      FFPvertex_->evaluate(ZERO,5,qin[ihel1].particle()->CC(),
				   qin[ihel1],pout2[phel2]);
	    diag[4] = FFGvertex_->evaluate(scale,inters6,inters4,gout[ghel]);
	    // eighth diagram
	    diag[5] = FFPvertex_->evaluate(ZERO,inters6,inters5,pout1[phel1]);
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // remove some coupling factors
  sum /= (8.*Constants::pi*generator()->standardModel()->alphaS(scale));
  // compute the two dipole terms
  double x = 1.-(pGqqbar_[order]*pQqqbar_[order]+pGqqbar_[order]*pQbarqqbar_[order])/
    (pQqqbar_[order]*pQbarqqbar_[order]);
  Lorentz5Momentum Kt = pQqqbar_[order]+pQbarqqbar_[order]-pGqqbar_[order];
  Lorentz5Momentum pa[4],pb[4];
  // momenta for first emission
  if(partons_[0]->id()>0) {
    pa[0] = x*pQqqbar_[order];
    pa[1] = pQbarqqbar_[order]; 
  }
  else {
    pa[0] = x*pQbarqqbar_[order];
    pa[1] = pQqqbar_[order]; 
  }
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=0;ix<2;++ix) {
    pa[ix+2] = pGammaqqbar_[order][ix]
      -2.*Ksum*(Ksum*pGammaqqbar_[order][ix])/Ksum2
      +2*K*(Kt*pGammaqqbar_[order][ix])/K2;
  }
  // momenta for second emission
  if(partons_[0]->id()>0) {
    pb[0] = pQqqbar_[order];
    pb[1] = x*pQbarqqbar_[order]; 
  }
  else {
    pb[0] = pQbarqqbar_[order];
    pb[1] = x*pQqqbar_[order]; 
  }
  K = pb[0]+pb[1];
  Ksum = K+Kt;
  K2 = K.m2();
  Ksum2 = Ksum.m2();
  for(unsigned int ix=0;ix<2;++ix) {
    pb[ix+2] = pGammaqqbar_[order][ix]
      -2.*Ksum*(Ksum*pGammaqqbar_[order][ix])/Ksum2
      +2*K*(Kt*pGammaqqbar_[order][ix])/K2;
  }
  // first LO matrix element
  Energy2 s,t,u;
  s = (pa[0]+pa[1]).m2();
  t = (pa[0]-pa[2]).m2();
  u = (pa[0]-pa[3]).m2();
  double lo1 = 8.*sqr(4.*Constants::pi*generator()
		      ->standardModel()->alphaEM(ZERO))*(t/u+u/t)*
    pow(double(partons_[0]->iCharge())/3.,4);
  // second LO matrix element
  s = (pb[0]+pb[1]).m2();
  t = (pb[0]-pb[2]).m2();
  u = (pb[0]-pb[3]).m2();
  double lo2 = 8.*sqr(4.*Constants::pi*generator()
		      ->standardModel()->alphaEM(ZERO))*(t/u+u/t)*
    pow(double(partons_[0]->iCharge())/3.,4);
  // first dipole
  InvEnergy2 D1 = 0.5/(pQqqbar_[order]*pGqqbar_[order])/x*(2./(1.-x)-(1.+x));
  // second dipole
  InvEnergy2 D2 = 0.5/(pQbarqqbar_[order]*pGqqbar_[order])/x*(2./(1.-x)-(1.+x));
  if(partons_[0]->id()<0) {
    swap(D1 ,D2 );
    swap(lo1,lo2);
  }
  // result
  double me;
  if(order==0) {
    //me = s_*abs(D1);
    me = s_*abs(D1)/(abs(D1)*lo1+abs(D2)*lo2)*UnitRemoval::InvE2*sum;
  }
  else {
    //me = s_*abs(D2);
    me = s_*abs(D2)/(abs(D1)*lo1+abs(D2)*lo2)*UnitRemoval::InvE2*sum;
  }
  return 4./3.*me;
}

void GammaGammaHardGenerator::generateQGQ() {
  // PDFs for the LO cross section
  double pdf[4];
  pdf[0]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],sqr(systemMass_),x_[0])/x_[0];
  pdf[1]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],sqr(systemMass_),x_[1])/x_[1];
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0;
  // storage of pt and rapidity of the quark
  Energy pT = 0.5*rs_;
  double yj,phi,x1,x2;
  // prefactor for veto algorithm
  double c = alphaS_->overestimateValue()/Constants::twopi*
    qgqFactor_*(maxyj-minyj);
  if(power_!=1.) c /= power_-1.;
  double wgt(0.);
  do {
    // generate pT
    if(power_!=1.) {
      pT = rs_*pow(pow(double(pT/rs_),1.-power_)
		   -log(UseRandom::rnd())/c,1./(1.-power_));
    }
    else {
      pT *= pow(UseRandom::rnd(),1./c);
    }
    // generate rapidity of the jet
    yj = UseRandom::rnd()*(maxyj-minyj)+ minyj;
    // generate phi
    phi = UseRandom::rnd()*Constants::twopi;
    // the CS variables
    double z,vt;
    if(partons_[0]->id()<0) {
      z  = (1.-pT/rs_/x_[1]*exp(-yj))/(1.+pT/rs_*exp( yj)/x_[0]);
      vt = pT/rs_/x_[1]*exp(-yj);
      x1 = x_[0]/z;
      x2 = x_[1];
      if(z<x_[0]||z>1.||vt>1.-z||vt<0) continue;
    }
    else {
      z = (1.-pT/rs_/x_[0]*exp( yj))/(1.+pT/rs_*exp(-yj)/x_[1]);
      vt = pT/rs_/x_[0]*exp( yj);
      x1 = x_[0];
      x2 = x_[1]/z;
      if(z<x_[1]||z>1.||vt>1.-z||vt<0) continue;
    }
    // momentum of the outgoing quark
    pQoutqg_ = Lorentz5Momentum (pT *cos(phi),pT *sin(phi),
				 pT *sinh(yj),pT *cosh(yj),ZERO);
    // momenta of the incoming quarks
    if(partons_[0]->id()>0) {
      pQinqg_ = Lorentz5Momentum (ZERO,ZERO,
				  x1*0.5*rs_,x1*0.5*rs_,ZERO);
      pGqg_   = Lorentz5Momentum (ZERO,ZERO,
				  -x2*0.5*rs_,x2*0.5*rs_,ZERO);
    }
    else {
      pGqg_   = Lorentz5Momentum (ZERO,ZERO,
				  x1*0.5*rs_,x1*0.5*rs_,ZERO);
      pQinqg_ = Lorentz5Momentum (ZERO,ZERO,
				  -x2*0.5*rs_,x2*0.5*rs_,ZERO);
    }
    // momenta of the photons
    Lorentz5Momentum Kt = loMomenta_[0]+loMomenta_[1];
    Lorentz5Momentum K  = pGqg_+pQinqg_-pQoutqg_;
    Lorentz5Momentum Ksum = K+Kt;
    Energy2 K2    = K   .m2();
    Energy2 Ksum2 = Ksum.m2();
    for(unsigned int iz=0;iz<2;++iz) { 
      pGammaqg_[iz] = loMomenta_[2+iz]
	-2.*(loMomenta_[2+iz]*Ksum)/Ksum2*Ksum+2.*(Kt*loMomenta_[2+iz])/K2*K;
    }
    // enforce LO photon cuts to avoid QED singularities
    if(pGammaqg_[0].perp() < minPhotonpT_  ||
       pGammaqg_[1].perp() < minPhotonpT_  ||
       pGammaqg_[0].eta()  < minPhotonEta_ ||
       pGammaqg_[1].eta()  < minPhotonEta_ ||
       pGammaqg_[0].eta()  > maxPhotonEta_ ||
       pGammaqg_[1].eta()  > maxPhotonEta_ ) {
      wgt=0.;
      continue;
    }
    // now for the weight
    // phase space/over estimate
    wgt  = 2./(1.-vt)*pow(pT/rs_,power_+1)/qgqFactor_;
    // pdf bit
    Energy2 scale = sqr(systemMass_)+sqr(pT);
    if(partons_[0]->id()>0) {
      pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x1)/x1;
      pdf[3]=beams_[1]->pdf()->xfx(beams_[1],gluon_     ,scale,x2)/x2;
    }
    else {
      pdf[2]=beams_[0]->pdf()->xfx(beams_[0],gluon_     ,scale,x1)/x1;
      pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,x2)/x2;
    }
    if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
      wgt=0.;
      continue;
    }
    wgt *= pdf[2]*pdf[3]/pdf[0]/pdf[1];
    // alpha S
    wgt *= alphaS_->ratio(sqr(pT));
    // final bit
    wgt *= QGQRatio();
    if(wgt>1.) {
      cerr << "testing weight QGQ " << wgt << " " << pT/GeV << " " << yj << "\n";
      cerr << "testing photon momenta " << pGammaqg_[0]/GeV << "\n";
      cerr << "testing photon momenta " << pGammaqg_[1]/GeV << "\n";
    }
  }
  while(UseRandom::rnd()>wgt&&pT>pTmin_);
  if(pT<pTmin_) pT=-GeV;
  pTqg_ = pT;
}

double GammaGammaHardGenerator::QGQRatio() {
  using namespace ThePEG::Helicity;
  double sum(0.);
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qout;
  vector<VectorWaveFunction> pout1,pout2,gin;
  SpinorWaveFunction q_in;
  SpinorBarWaveFunction q_out;
  if(partons_[0]->id()>0) {
    q_in  = SpinorWaveFunction   (pQinqg_ ,partons_[0]      ,incoming);
    q_out = SpinorBarWaveFunction(pQoutqg_,partons_[1]->CC(),outgoing);
  }
  else {
    q_in  = SpinorWaveFunction   (pQinqg_ ,partons_[1]      ,incoming);
    q_out = SpinorBarWaveFunction(pQoutqg_,partons_[0]->CC(),outgoing);
  }
  VectorWaveFunction p_out1(pGammaqg_[0],photon_,outgoing);
  VectorWaveFunction p_out2(pGammaqg_[1],photon_,outgoing);
  VectorWaveFunction g_in  (pGqg_       ,gluon_ ,incoming);
  for(unsigned int ix=0;ix<2;++ix) {
    q_in .reset(ix);         qin  .push_back(q_in  );
    q_out.reset(ix);         qout .push_back(q_out );
    g_in .reset(2*ix);       gin  .push_back(g_in  );
    p_out1.reset(2*ix);      pout1.push_back(p_out1);
    p_out2.reset(2*ix);      pout2.push_back(p_out2);
  }
  vector<Complex> diag(6);
  Energy2 scale = sqr(systemMass_);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int phel1=0;phel1<2;++phel1) {
	for(unsigned int phel2=0;phel2<2;++phel2) {
	  for(unsigned int ohel=0;ohel<2;++ohel) {
 	    // first diagram
	    SpinorWaveFunction inters1 = 
	      FFPvertex_->evaluate(ZERO,5,qin[ihel1].particle()->CC(),
				   qin[ihel1],pout1[phel1]);
	    SpinorBarWaveFunction inters2 = 
	      FFPvertex_->evaluate(ZERO,5,qout[ohel].particle(),
				   qout[ohel],pout2[phel2]);
	    diag[0] = FFGvertex_->evaluate(scale,inters1,inters2,gin[ihel2]);
	    // second diagram
	    SpinorWaveFunction inters3 = 
	      FFGvertex_->evaluate(scale,5,qin[ihel1].particle()->CC(),
				   qin[ihel1],gin[ihel2]);
	    SpinorBarWaveFunction inters4 = 
	      FFPvertex_->evaluate(ZERO,5,qout[ohel].particle(),
				   qout[ohel],pout1[phel1]);
	    diag[1] = FFPvertex_->evaluate(ZERO,inters3,inters4,pout2[phel2]);
	    // fourth diagram
	    diag[2] = FFPvertex_->evaluate(ZERO,inters3,inters2,pout1[phel1]);
	    // fifth diagram
	    SpinorBarWaveFunction inters5 = 
	      FFGvertex_->evaluate(scale,5,qout[ohel].particle(),
				   qout[ohel],gin[ihel2]);
	    diag[3] = 
	      FFPvertex_->evaluate(ZERO,inters1,inters5,pout2[phel2]);
	    // sixth diagram
	    SpinorWaveFunction inters6 = 
	      FFPvertex_->evaluate(ZERO,5,qin[ihel1].particle()->CC(),
				   qin[ihel1],pout2[phel2]);
	    diag[4] = FFGvertex_->evaluate(scale,inters6,inters4,gin[ihel2]);
	    // eighth diagram
	    diag[5] = FFPvertex_->evaluate(ZERO,inters6,inters5,pout1[phel1]);
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // remove some coupling factors
  sum /= (8.*Constants::pi*generator()->standardModel()->alphaS(scale));
  // compute the two dipole terms
  double x = 1.-(pGqg_*pQoutqg_+pQinqg_*pQoutqg_)/(pGqg_*pQinqg_);
  Lorentz5Momentum Kt = pQinqg_+pGqg_-pQoutqg_;
  Lorentz5Momentum pa[4],pb[4],pc[4];
  // momenta for IS emission
  if(partons_[0]->id()>0) {
    pa[0] = pQinqg_;
    pa[1] = x*pGqg_; 
  }
  else {
    pa[0] = x*pGqg_;
    pa[1] = pQinqg_; 
  }
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=0;ix<2;++ix) {
    pa[ix+2] = pGammaqg_[ix]
      -2.*Ksum*(Ksum*pGammaqg_[ix])/Ksum2
      +2*K*(Kt*pGammaqg_[ix])/K2;
  }
  // momenta for first FS emission
  double xFS[2],zFS[2];
  for(unsigned int ix=0;ix<2;++ix) {
    xFS[ix] = 1.-(pGammaqg_[ix]*pQoutqg_)/((pGammaqg_[ix]+pQoutqg_)*pQinqg_);
    zFS[ix] = (pQoutqg_*pQinqg_)/((pGammaqg_[ix]+pQoutqg_)*pQinqg_);
  }
  // first set of momenta
  pb[0] = xFS[0]*pQinqg_;
  pb[1] = pGqg_;
  pb[2] = pGammaqg_[1];
  pb[3] = pGammaqg_[0]+pQoutqg_-(1.-xFS[0])*pQinqg_;
  // second set of momenta
  pc[0] = xFS[1]*pQinqg_;
  pc[1] = pGqg_;
  pc[2] = pGammaqg_[0];
  pc[3] = pGammaqg_[1]+pQoutqg_-(1.-xFS[1])*pQinqg_;
  // first LO matrix element
  Energy2 s,t,u;
  s = (pa[0]+pa[1]).m2();
  t = (pa[0]-pa[2]).m2();
  u = (pa[0]-pa[3]).m2();
  double coupling = 
    sqr(4.*Constants::pi*generator()->standardModel()->alphaEM(ZERO))*
    pow(double(partons_[0]->iCharge())/3.,4);
  double lo1 = 8.*coupling*(t/u+u/t);
  // second LO matrix element
  s = (pb[0]+pb[1]).m2();
  t = (pb[0]-pb[2]).m2();
  u = (pb[0]-pb[3]).m2();
  double lo2 = -8./s/t*(s*s+t*t+2.*u*(s+t+u))*coupling;
  // third  LO matrix element
  s = (pc[0]+pc[1]).m2();
  t = (pc[0]-pc[2]).m2();
  u = (pc[0]-pc[3]).m2();
  double lo3 = -8./s/t*(s*s+t*t+2.*u*(s+t+u))*coupling;
  // first dipole
  InvEnergy2 D1 =  0.5/(pGqg_*pQoutqg_)/x*(1.-2.*x*(1.-x));
  // second dipole
  InvEnergy2 D2 =  0.5/(pQoutqg_*pGammaqg_[0])*(2./(2.-zFS[0]-xFS[0])-1.-zFS[0]);
  // third  dipole
  InvEnergy2 D3 =  0.5/(pQoutqg_*pGammaqg_[1])*(2./(2.-zFS[1]-xFS[1])-1.-zFS[1]);
  // result
  return s_*abs(D1)/(abs(D1)*lo1+abs(D2)*lo2+abs(D3)*lo3)*UnitRemoval::InvE2*sum;
}

void GammaGammaHardGenerator::generateQbarGQbar() {
  // PDFs for the LO cross section
  double pdf[4];
  pdf[0]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],sqr(systemMass_),x_[0])/x_[0];
  pdf[1]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],sqr(systemMass_),x_[1])/x_[1];
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0;
  // storage of pt and rapidity of the quark
  Energy pT = 0.5*rs_;
  double yj,phi,x1,x2;
  // prefactor for veto algorithm
  double c = alphaS_->overestimateValue()/Constants::twopi*
    qbargqbarFactor_*(maxyj-minyj);
  if(power_!=1.) c /= power_-1.;
  double wgt(0.);
  do {
    // generate pT
    if(power_!=1.) {
      pT = rs_*pow(pow(double(pT/rs_),1.-power_)
		   -log(UseRandom::rnd())/c,1./(1.-power_));
    }
    else {
      pT *= pow(UseRandom::rnd(),1./c);
    }
    // generate rapidity of the jet
    yj = UseRandom::rnd()*(maxyj-minyj)+ minyj;
    // generate phi
    phi = UseRandom::rnd()*Constants::twopi;
    // the CS variables
    double z,vt;
    if(partons_[0]->id()>0) {
      z  = (1.-pT/rs_/x_[1]*exp(-yj))/(1.+pT/rs_*exp( yj)/x_[0]);
      vt = pT/rs_/x_[1]*exp(-yj);
      x1 = x_[0]/z;
      x2 = x_[1];
      if(z<x_[0]||z>1.||vt>1.-z||vt<0) continue;
    }
    else {
      z = (1.-pT/rs_/x_[0]*exp( yj))/(1.+pT/rs_*exp(-yj)/x_[1]);
      vt = pT/rs_/x_[0]*exp( yj);
      x1 = x_[0];
      x2 = x_[1]/z;
      if(z<x_[1]||z>1.||vt>1.-z||vt<0) continue;
    }
    // momentum of the outgoing quark
    pQoutgqbar_ = Lorentz5Momentum (pT *cos(phi),pT *sin(phi),
				    pT *sinh(yj),pT *cosh(yj),ZERO);
    // momenta of the incoming quarks
    if(partons_[0]->id()<0) {
      pQingqbar_ = Lorentz5Momentum (ZERO,ZERO,
				     x1*0.5*rs_,x1*0.5*rs_,ZERO);
      pGgqbar_   = Lorentz5Momentum (ZERO,ZERO,
				     -x2*0.5*rs_,x2*0.5*rs_,ZERO);
    }
    else {
      pGgqbar_   = Lorentz5Momentum (ZERO,ZERO,
				     x1*0.5*rs_,x1*0.5*rs_,ZERO);
      pQingqbar_ = Lorentz5Momentum (ZERO,ZERO,
				     -x2*0.5*rs_,x2*0.5*rs_,ZERO);
    }
    // momenta of the photons
    Lorentz5Momentum Kt = loMomenta_[0]+loMomenta_[1];
    Lorentz5Momentum K  = pGgqbar_+pQingqbar_-pQoutgqbar_;
    Lorentz5Momentum Ksum = K+Kt;
    Energy2 K2    = K   .m2();
    Energy2 Ksum2 = Ksum.m2();
    for(unsigned int iz=0;iz<2;++iz) { 
      pGammagqbar_[iz] = loMomenta_[2+iz]
	-2.*(loMomenta_[2+iz]*Ksum)/Ksum2*Ksum+2.*(Kt*loMomenta_[2+iz])/K2*K;
    }
    // enforce LO photon cuts to avoid QED singularities
    if(pGammagqbar_[0].perp() < minPhotonpT_  ||
       pGammagqbar_[1].perp() < minPhotonpT_  ||
       pGammagqbar_[0].eta()  < minPhotonEta_ ||
       pGammagqbar_[1].eta()  < minPhotonEta_ ||
       pGammagqbar_[0].eta()  > maxPhotonEta_ ||
       pGammagqbar_[1].eta()  > maxPhotonEta_ ) {
      wgt=0.;
      continue;
    }
    // now for the weight
    // phase space/over estimate
    wgt  = 2./(1.-vt)*pow(pT/rs_,power_+1)/qbargqbarFactor_;
    // pdf bit
    Energy2 scale = sqr(systemMass_)+sqr(pT);
    if(partons_[0]->id()<0) {
      pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x1)/x1;
      pdf[3]=beams_[1]->pdf()->xfx(beams_[1],gluon_     ,scale,x2)/x2;
    }
    else {
      pdf[2]=beams_[0]->pdf()->xfx(beams_[0],gluon_     ,scale,x1)/x1;
      pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,x2)/x2;
    }
    if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
      wgt=0.;
      continue;
    }
    wgt *= pdf[2]*pdf[3]/pdf[0]/pdf[1];
    // alpha S
    wgt *= alphaS_->ratio(sqr(pT));
    // final bit
    wgt *= QbarGQbarRatio();
    if(wgt>1.) {
      cerr << "testing weight GQBARQ " << wgt << " " << pT/GeV << " " << yj << "\n";
      cerr << "testing photon momenta " << pGammagqbar_[0]/GeV << "\n";
      cerr << "testing photon momenta " << pGammagqbar_[1]/GeV << "\n";
    }
  }
  while(UseRandom::rnd()>wgt&&pT>pTmin_);
  if(pT<pTmin_) pT=-GeV;
  pTgqbar_ = pT;
}

double GammaGammaHardGenerator::QbarGQbarRatio() {
  using namespace ThePEG::Helicity;
  double sum(0.);
  vector<SpinorBarWaveFunction> qin;
  vector<SpinorWaveFunction> qout;
  vector<VectorWaveFunction> pout1,pout2,gin;
  SpinorBarWaveFunction q_in;
  SpinorWaveFunction q_out;
  if(partons_[0]->id()<0) {
    q_in  = SpinorBarWaveFunction(pQingqbar_ ,partons_[0]      ,incoming);
    q_out = SpinorWaveFunction   (pQoutgqbar_,partons_[1]->CC(),outgoing);
  }
  else {
    q_in  = SpinorBarWaveFunction(pQingqbar_ ,partons_[1]      ,incoming);
    q_out = SpinorWaveFunction   (pQoutgqbar_,partons_[0]->CC(),outgoing);
  }
  VectorWaveFunction p_out1(pGammagqbar_[0],photon_,outgoing);
  VectorWaveFunction p_out2(pGammagqbar_[1],photon_,outgoing);
  VectorWaveFunction g_in  (pGgqbar_       ,gluon_ ,incoming);
  for(unsigned int ix=0;ix<2;++ix) {
    q_in .reset(ix);         qin  .push_back(q_in  );
    q_out.reset(ix);         qout .push_back(q_out );
    g_in .reset(2*ix);       gin  .push_back(g_in  );
    p_out1.reset(2*ix);      pout1.push_back(p_out1);
    p_out2.reset(2*ix);      pout2.push_back(p_out2);
  }
  vector<Complex> diag(6);
  Energy2 scale = sqr(systemMass_);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int phel1=0;phel1<2;++phel1) {
	for(unsigned int phel2=0;phel2<2;++phel2) {
	  for(unsigned int ohel=0;ohel<2;++ohel) {
  	    // first diagram
	    SpinorBarWaveFunction inters1 = 
	      FFPvertex_->evaluate(ZERO,5,qin[ihel1].particle()->CC(),
 				   qin[ihel1],pout1[phel1]);
	    SpinorWaveFunction inters2 = 
	      FFPvertex_->evaluate(ZERO,5,qout[ohel].particle(),
				   qout[ohel],pout2[phel2]);
	    diag[0] = FFGvertex_->evaluate(scale,inters2,inters1,gin[ihel2]);
	    // second diagram
	    SpinorBarWaveFunction inters3 = 
	      FFGvertex_->evaluate(scale,5,qin[ihel1].particle()->CC(),
				   qin[ihel1],gin[ihel2]);
	    SpinorWaveFunction inters4 = 
	      FFPvertex_->evaluate(ZERO,5,qout[ohel].particle(),
				   qout[ohel],pout1[phel1]);
	    diag[1] = FFPvertex_->evaluate(ZERO,inters4,inters3,pout2[phel2]);
	    // fourth diagram
	    diag[2] = FFPvertex_->evaluate(ZERO,inters2,inters3,pout1[phel1]);
	    // fifth diagram
	    SpinorWaveFunction inters5 = 
	      FFGvertex_->evaluate(scale,5,qout[ohel].particle(),
				   qout[ohel],gin[ihel2]);
	    diag[3] = 
	      FFPvertex_->evaluate(ZERO,inters5,inters1,pout2[phel2]);
	    // sixth diagram
	    SpinorBarWaveFunction inters6 = 
	      FFPvertex_->evaluate(ZERO,5,qin[ihel1].particle()->CC(),
				   qin[ihel1],pout2[phel2]);
	    diag[4] = FFGvertex_->evaluate(scale,inters4,inters6,gin[ihel2]);
	    // eighth diagram
	    diag[5] = FFPvertex_->evaluate(ZERO,inters5,inters6,pout1[phel1]);
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // remove some coupling factors
  sum /= (8.*Constants::pi*generator()->standardModel()->alphaS(scale));
  // compute the two dipole terms
  double x = 1.-(pGgqbar_*pQoutgqbar_+pQingqbar_*pQoutgqbar_)/(pGgqbar_*pQingqbar_);
  Lorentz5Momentum Kt = pQingqbar_+pGgqbar_-pQoutgqbar_;
  Lorentz5Momentum pa[4],pb[4],pc[4];
  // momenta for IS emission
  if(partons_[0]->id()<0) {
    pa[0] = pQingqbar_;
    pa[1] = x*pGgqbar_; 
  }
  else {
    pa[0] = x*pGgqbar_;
    pa[1] = pQingqbar_; 
  }
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=0;ix<2;++ix) {
    pa[ix+2] = pGammagqbar_[ix]
      -2.*Ksum*(Ksum*pGammagqbar_[ix])/Ksum2
      +2*K*(Kt*pGammagqbar_[ix])/K2;
  }
  // momenta for first FS emission
  double xFS[2],zFS[2];
  for(unsigned int ix=0;ix<2;++ix) {
    xFS[ix] = 1.-(pGammagqbar_[ix]*pQoutgqbar_)/
      ((pGammagqbar_[ix]+pQoutgqbar_)*pQingqbar_);
    zFS[ix] = (pQoutgqbar_*pQingqbar_)/((pGammagqbar_[ix]+pQoutgqbar_)*pQingqbar_);
  }
  // first set of momenta
  pb[0] = pGgqbar_;
  pb[1] = xFS[0]*pQingqbar_;
  pb[2] = pGammagqbar_[1];
  pb[3] = pGammagqbar_[0]+pQoutgqbar_-(1.-xFS[0])*pQingqbar_;
  // second set of momenta
  pc[0] = pGgqbar_;
  pc[1] = xFS[1]*pQingqbar_;
  pc[2] = pGammagqbar_[0];
  pc[3] = pGammagqbar_[1]+pQoutgqbar_-(1.-xFS[1])*pQingqbar_;
  // first LO matrix element
  Energy2 s,t,u;
  s = (pa[0]+pa[1]).m2();
  t = (pa[0]-pa[2]).m2();
  u = (pa[0]-pa[3]).m2();
  double coupling = 
    sqr(4.*Constants::pi*generator()->standardModel()->alphaEM(ZERO))*
    pow(double(partons_[0]->iCharge())/3.,4);
  double lo1 = 8.*coupling*(t/u+u/t);
  // second LO matrix element
  s = (pb[0]+pb[1]).m2();
  u = (pb[0]-pb[2]).m2();
  t = (pb[0]-pb[3]).m2();
  double lo2 = -8./s/t*(s*s+t*t+2.*u*(s+t+u))*coupling;
  // third  LO matrix element
  s = (pc[0]+pc[1]).m2();
  u = (pc[0]-pc[2]).m2();
  t = (pc[0]-pc[3]).m2();
  double lo3 = -8./s/t*(s*s+t*t+2.*u*(s+t+u))*coupling;
  // first dipole
  InvEnergy2 D1 =  0.5/(pGgqbar_*pQoutgqbar_)/x*(1.-2.*x*(1.-x));
  // second dipole
  InvEnergy2 D2 =  0.5/(pQoutgqbar_*pGammagqbar_[0])*(2./(2.-zFS[0]-xFS[0])-1.-zFS[0]);
  // third  dipole
  InvEnergy2 D3 =  0.5/(pQoutgqbar_*pGammagqbar_[1])*(2./(2.-zFS[1]-xFS[1])-1.-zFS[1]);
  // result
  return s_*abs(D1)/(abs(D1)*lo1+abs(D2)*lo2+abs(D3)*lo3)*UnitRemoval::InvE2*sum;
}
