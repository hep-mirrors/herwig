// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VGammaHardGenerator class.
//

#include "VGammaHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <numeric>


using namespace Herwig;

VGammaHardGenerator::VGammaHardGenerator() : pTmin_(2.*GeV), qqgFactor_(1.), power_(2.) {}

IBPtr VGammaHardGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr VGammaHardGenerator::fullclone() const {
  return new_ptr(*this);
}

void VGammaHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << alphaS_ << gluon_ << ounit(pTmin_,GeV)
     << FFPvertex_ << FFWvertex_ << FFZvertex_ << WWWvertex_ << FFGvertex_;
}

void VGammaHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> alphaS_ >> gluon_ >> iunit(pTmin_,GeV)
     >> FFPvertex_ >> FFWvertex_ >> FFZvertex_ >> WWWvertex_ >> FFGvertex_;
}

ClassDescription<VGammaHardGenerator> VGammaHardGenerator::initVGammaHardGenerator;
// Definition of the static class description member.

void VGammaHardGenerator::Init() {

  static ClassDocumentation<VGammaHardGenerator> documentation
    ("The VGammaHardGenerator class implements the generation of the "
     "hard QCD radiation in electroweak vector boson+photon processes");

  static Reference<VGammaHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &VGammaHardGenerator::alphaS_, false, false, true, false, false);

  static Parameter<VGammaHardGenerator, Energy> interfacePtMin
    ("minPt",
     "The pt cut on hardest emision generation"
     "2*(1-Beta)*exp(-sqr(intrinsicpT/RMS))/sqr(RMS)",
     &VGammaHardGenerator::pTmin_, GeV, 2.*GeV, ZERO, 100000.0*GeV,
     false, false, Interface::limited);
}

void VGammaHardGenerator::doinit() {
  HardestEmissionGenerator::doinit();
  gluon_ = getParticleData(ParticleID::g);
  // get a pointer to the standard model object in the run
  static const tcHwSMPtr hwsm
    = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if (!hwsm) throw InitException() << "hwsm pointer is null in"
				   << " MEPP2VGamma::doinit()"
				   << Exception::abortnow;
  // get pointers to all required Vertex objects
  FFZvertex_ = hwsm->vertexFFZ();
  FFPvertex_ = hwsm->vertexFFP();
  WWWvertex_ = hwsm->vertexWWW();
  FFWvertex_ = hwsm->vertexFFW();
  FFGvertex_ = hwsm->vertexFFG();
}

bool VGammaHardGenerator::canHandle(ShowerTreePtr tree) {
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
  ix=0;  
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    part[ix]=cjt->first->progenitor();
    ++ix;
  }
  // put photon first
  if(part[0]->id()!=ParticleID::gamma) swap(part[0],part[1]);
  // check photon
  if(part[0]->id()!=ParticleID::gamma) return false;
  // check gauge boson
  if(abs(part[1]->id())==ParticleID::Wplus ||
     part[1]->id() ==ParticleID::Z0) return true;
  else return false;
}

HardTreePtr VGammaHardGenerator::generateHardest(ShowerTreePtr tree) {
  // get the particles to be showered
  beams_.clear();
  partons_.clear();
  // find the incoming particles
  ShowerParticleVector incoming;
  quarkplus_ = true;
  vector<ShowerProgenitorPtr> particlesToShower;
  //progenitor particles are produced in z direction.
  for( map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	 cit = tree->incomingLines().begin(); 
       cit != tree->incomingLines().end(); ++cit ) {
    incoming.push_back( cit->first->progenitor() );
    beams_.push_back( cit->first->beam() );
    partons_.push_back( cit->first->progenitor()->dataPtr() );
    // check that quark is along +ve z direction
    if(cit->first->progenitor()->id() > 0 && 
       cit->first->progenitor()->momentum().z() < ZERO ) 
      quarkplus_ = false;
    particlesToShower.push_back( cit->first );
  }
  // find the outgoing particles
  tShowerParticlePtr boson,photon;
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
 	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    if(cjt->first->progenitor()->id()==ParticleID::gamma)
      photon = cjt->first->progenitor();
    else if(abs(cjt->first->progenitor()->id())==ParticleID::Wplus||
	    cjt->first->progenitor()->id()==ParticleID::Z0)
      boson  =  cjt->first->progenitor();
    particlesToShower.push_back(cjt->first);
  }
  assert(boson&&photon);
  // we are assuming quark first, swap order to ensure this
  // if antiquark first
  if(partons_[0]->id()<partons_[1]->id()) {
    swap(partons_[0],partons_[1]);
    swap(beams_[0],beams_[1]);
  }
  // Born variables
  // Rapidity of the photon
  photonRapidity_ = photon->momentum().rapidity();
  // Rapidity of the gauge boson
  bosonRapidity_  =  boson->momentum().rapidity();
  // pT of the photon
  photonpT_       = photon->momentum().perp();
  // Azimuth of the photon
  photonAzimuth_  = photon->momentum().phi();
  // gauge boson mass
  bosonMass_      = boson->mass();
  // mass of the boson/photon system
  systemMass_     = (photon->momentum()+boson->momentum()).m();
  // direction flip if needed
  if(!quarkplus_) {
    photonRapidity_ *= -1.;
    bosonRapidity_  *= -1.;
    photonAzimuth_  *= -1.;
  }
  // ParticleData objects
  photon_ = photon->dataPtr();
  boson_  = boson->dataPtr();
  // CMS energy of the hadron collision
  s_  = generator()->currentEvent()->primaryCollision()->m2();
  rs_ = sqrt(s_);
  // Borm momentum fractions
  double systemRapidity = (photon->momentum()+boson->momentum()).rapidity();
  if(!quarkplus_) systemRapidity *= -1.;
  x_[0] = systemMass_*exp( systemRapidity)/sqrt(s_);          
  x_[1] = systemMass_*exp(-systemRapidity)/sqrt(s_);
  // generate the new configuration
  pTqqbar_ = -GeV;
  pTqg_    = -GeV;
  pTgqbar_ = -GeV;
  // generate a q qbar -> V g cconfiguration
  generateQQbarG();
  // return if no emission
  if(pTqqbar_<ZERO&&pTqg_<ZERO&&pTgqbar_<ZERO) {
    for(unsigned int ix=0;ix<particlesToShower.size();++ix)
      particlesToShower[ix]->maximumpT(pTmin_);
    return HardTreePtr();
  }
  // select the type of emission
  int emissionType=0;
  if      (pTqqbar_>pTqg_    && pTqqbar_>pTgqbar_ ) emissionType = 1;
  else if (pTqg_>pTqqbar_    && pTqg_>pTgqbar_    ) emissionType = 2;
  else if (pTgqbar_>pTqqbar_ && pTgqbar_>pTqg_    ) emissionType = 3;
  int iemit;
  ShowerParticleVector newparticles;
  Lorentz5Momentum pboson,pphoton;
  Energy pTEmit;
  // q qbar -> V gamma g
  if(emissionType==1) {
    newparticles.push_back(new_ptr(ShowerParticle(partons_[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(partons_[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gluon_            , true)));
    newparticles[0]->set5Momentum(pQqqbar_);
    newparticles[1]->set5Momentum(pQbarqqbar_);
    newparticles[2]->set5Momentum(pGqqbar_);
    pboson  = pVqqbar_;
    pphoton = pGammaqqbar_; 
    iemit = pQqqbar_.z()/pGqqbar_.rapidity()>ZERO ? 0 : 1;
    pTEmit = pTqqbar_;
  }
  // q g    -> q V gamma
  else if(emissionType==2) {
    iemit=1;
    newparticles.push_back(new_ptr(ShowerParticle(partons_[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gluon_           ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(partons_[1]->CC(), true)));
    newparticles[0]->set5Momentum(pQinqg_);
    newparticles[1]->set5Momentum(pGqg_);
    newparticles[2]->set5Momentum(pQoutqg_);
    pboson  = pVqg_;
    pphoton = pGammaqg_;
    pTEmit = pTqg_;
  }
  // g qbar -> qbar V gamma
  else {
    iemit=0;
    newparticles.push_back(new_ptr(ShowerParticle(gluon_           ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(partons_[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(partons_[0]->CC(), true)));
    newparticles[0]->set5Momentum(pGgqbar_);
    newparticles[1]->set5Momentum(pQingqbar_);
    newparticles[2]->set5Momentum(pQoutgqbar_);
    pboson  = pVgqbar_;
    pphoton = pGammagqbar_;
    pTEmit  = pTgqbar_;
  }
  // create the photon
  newparticles.push_back(new_ptr(ShowerParticle(photon->dataPtr(),true)));
  newparticles.back()->set5Momentum(pphoton);
  // create the boson
  newparticles.push_back(new_ptr(ShowerParticle(boson ->dataPtr(),true)));
  newparticles.back()->set5Momentum(pboson);
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
				 << "VGammaHardGenerator::generateHardest()" 
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
  // outgoing boson
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
  for(set<HardBranchingPtr>::const_iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      if((((**cit).status()==HardBranching::Incoming && 
	   !particlesToShower[ix]->progenitor()->isFinalState())||
	  ((**cit).status()==HardBranching::Outgoing &&
	   particlesToShower[ix]->progenitor()->isFinalState()))&&
	 particlesToShower[ix]->progenitor()->id()==(**cit).branchingParticle()->id()) {
	particlesToShower[ix]->progenitor()->set5Momentum((**cit).showerMomentum());
      }
    }
  }
  return hardTree;
}

void VGammaHardGenerator::generateQQbarG() {
  Energy2 sHat;
  Energy rsHat;
  // storage of pt and rapidity of the gluon
  Energy pT = 0.5*rs_,pTV,mTV;
  double yj,phi,x1,x2,phiV,wgt;
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0;
  // PDFs for the LO cross section
  double pdf[4];
  pdf[0]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],sqr(systemMass_),x_[0]);
  pdf[1]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],sqr(systemMass_),x_[1]);
  // prefactor for veto algorithm
  double c = alphaS_->overestimateValue()/Constants::twopi*
    qqgFactor_*(maxyj-minyj)/(power_-1.);
  do {
    // generate pT
    pT = rs_*pow(pow(double(pT/rs_),1.-power_)-log(UseRandom::rnd())/c,1./(1.-power_));
    // generate rapidity of the jet
    yj = UseRandom::rnd()*(maxyj-minyj)+ minyj;
    // generate phi
    phi = UseRandom::rnd()*Constants::twopi;
    // pT of the W/Z
    pTV = sqrt(sqr(photonpT_*sin(photonAzimuth_)+pT*sin(phi))+
	       sqr(photonpT_*cos(photonAzimuth_)+pT*cos(phi)));
    mTV = sqrt(sqr(pTV)+sqr(bosonMass_));
    // azimuth of W/Z
    phiV = Constants::pi+atan2(photonpT_*sin(photonAzimuth_)+pT*sin(phi),
			       photonpT_*cos(photonAzimuth_)+pT*cos(phi));
    // calculate x_1 and x_2
    x1 = 0.5/rs_*(photonpT_*exp( photonRapidity_)+mTV*exp( bosonRapidity_)+pT*exp( yj));
    x2 = 0.5/rs_*(photonpT_*exp(-photonRapidity_)+mTV*exp(-bosonRapidity_)+pT*exp(-yj));
    // sHat
    sHat  = x1*x2*s_;
    rsHat = sqrt(sHat);
    // now for the weight
    // first kinematical factors
    wgt = 0.5*pow<4,1>(systemMass_)/sqr(sHat);
    // pdf bit
    Energy2 scale = sqr(systemMass_)+sqr(pT);
    pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x1);
    pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,x2);
    if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
      wgt=0.;
      continue;
    }
    wgt *= pdf[2]*pdf[3]/pdf[0]/pdf[1];
    // divide by most of over estimate
    wgt *= pow(pT/rs_,power_+1)/qqgFactor_;
    pGammaqqbar_ = Lorentz5Momentum (photonpT_*cos(photonAzimuth_),
				     photonpT_*sin(photonAzimuth_),
				     photonpT_*sinh(photonRapidity_),
				     photonpT_*cosh(photonRapidity_),ZERO);
    pVqqbar_     = Lorentz5Momentum (pTV*cos(phiV),pTV*sin(phiV),
				     mTV*sinh(bosonRapidity_),
				     mTV*cosh(bosonRapidity_),bosonMass_);
    pGqqbar_     = Lorentz5Momentum (pT *cos(phi ),pT *sin(phi ),
				     pT *sinh(       yj     ),
				     pT *cosh(      yj      ),ZERO);
    pQqqbar_     = Lorentz5Momentum (ZERO,ZERO, x1*rs_,x1*rs_,ZERO);
    pQbarqqbar_  = Lorentz5Momentum (ZERO,ZERO,-x2*rs_,x2*rs_,ZERO);
    // final bit
    wgt *= QQbarGratio();
  }
  while(UseRandom::rnd()>wgt&&pT>pTmin_);
  if(pT<pTmin_) pT=-GeV;
  pTqqbar_ = pT;
  if(!quarkplus_) {
    LorentzRotation trans;
    trans.rotateX(Constants::pi);
    pGammaqqbar_.transform(trans);
    pVqqbar_    .transform(trans);
    pGqqbar_    .transform(trans);
    pQqqbar_    .transform(trans);
    pQbarqqbar_ .transform(trans);
  }
  generator()->log() << "testing new \n"
		     << "Photon    = " << pGammaqqbar_/GeV << "\n"
		     << "Boson     = " << pVqqbar_/GeV << "\n"
		     << "Gluon     = " << pGqqbar_/GeV << "\n"
		     << "Quark     = " << pQqqbar_/GeV << "\n"
		     << "Antiquark = " << pQbarqqbar_/GeV << "\n"
		     << "sum " 
		     << (pGammaqqbar_+pVqqbar_+pGqqbar_-pQqqbar_-pQbarqqbar_)/GeV << "\n";
}

double VGammaHardGenerator::QQbarGratio() {
  using namespace ThePEG::Helicity;
  double sum(0.);
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qbarin;
  vector<VectorWaveFunction> wout,pout,gout;
  SpinorWaveFunction q_in(pQqqbar_,partons_[0],incoming);
  SpinorBarWaveFunction qbar_in(pQbarqqbar_,partons_[1],incoming);
  VectorWaveFunction v_out(pVqqbar_    ,boson_,outgoing);
  VectorWaveFunction p_out(pGammaqqbar_,photon_,outgoing);
  VectorWaveFunction g_out(pGqqbar_    ,gluon_,outgoing);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<2) {
      q_in.reset(ix);
      qin.push_back(q_in);
      qbar_in.reset(ix);
      qbarin.push_back(qbar_in);
      //g_out.reset(10);
      g_out.reset(2*ix);
      gout.push_back(g_out);
      //p_out.reset(10);
      p_out.reset(2*ix);
      pout.push_back(p_out);
    }
    v_out.reset(ix);
    wout.push_back(v_out);
  }
  Energy2 scale(sqr(bosonMass_));
  vector<Complex> diag(boson_->id()==ParticleID::Z0 ? 6 : 8, 0.);
  AbstractFFVVertexPtr vertex = 
    boson_->id()==ParticleID::Z0 ? FFZvertex_ : FFWvertex_;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int whel=0;whel<3;++whel) {
	for(unsigned int phel=0;phel<2;++phel) {
	  for(unsigned int ghel=0;ghel<2;++ghel) {
	    //
	    //  Diagrams which are the same for W/Z gamma
	    //
	    // first diagram
	    SpinorWaveFunction inters1 = 
	      FFPvertex_->evaluate(scale,5,partons_[0],qin[ihel1],pout[phel]);
	    SpinorBarWaveFunction inters2 = 
	      vertex->evaluate(scale,5,partons_[0]->CC(),qbarin[ihel2],wout[whel]);
	    diag[0] = FFGvertex_->evaluate(scale,inters1,inters2,gout[ghel]);
	    // second diagram
	    SpinorWaveFunction inters3 = 
	      FFGvertex_->evaluate(scale,5,partons_[0],qin[ihel1],gout[ghel]);
	    SpinorBarWaveFunction inters4 = 
	      FFPvertex_->evaluate(scale,5,partons_[1],qbarin[ihel2],pout[phel]);
	    diag[1] = vertex->evaluate(scale,inters3,inters4,wout[whel]);
	    // fourth diagram
	    diag[2] = FFPvertex_->evaluate(scale,inters3,inters2,pout[phel]);
	    // fifth diagram
	    SpinorBarWaveFunction inters5 = 
	      FFGvertex_->evaluate(scale,5,partons_[1],qbarin[ihel2],gout[ghel]);
	    diag[3] = 
	      vertex->evaluate(scale,inters1,inters5,wout[whel]);
	    // sixth diagram
	    SpinorWaveFunction inters6 = 
	      vertex->evaluate(scale,5,partons_[1]->CC(),qin[ihel1],wout[whel]);
	    diag[4] = FFGvertex_->evaluate(scale,inters6,inters4,gout[ghel]);
	    // eighth diagram
	    diag[5] = FFPvertex_->evaluate(scale,inters6,inters5,pout[phel]);
	    //
	    //  Diagrams only for W gamma
	    //
	    if(boson_->id()!=ParticleID::Z0) {
	      // third diagram
	      VectorWaveFunction interv = 
		WWWvertex_->evaluate(scale,3,boson_->CC(),pout[phel],wout[whel]);
	      diag[6] = vertex->evaluate(scale,inters3,qbarin[ihel2],interv);
	      // seventh diagram
	      diag[7] = vertex->evaluate(scale,qin[ihel1],inters5,interv);
	    }
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    cerr << "testing" 
		 << ihel1 << " " << ihel2 << " " 
		 << ghel  << " " << phel  << " " << whel << " " 
		 << dsum/sqrt(norm(diag[0])+norm(diag[1])+norm(diag[2])+norm(diag[3])+
			      norm(diag[4])+norm(diag[5])+norm(diag[6])+norm(diag[7]))<< "\n";
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // final spin and colour factors
  // spin = 1/4 colour = 4/9
  sum /= 9.;


  return 1000000.;
}
