// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVHardGenerator class.
//
#include <math.h>

#include "VVHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/Hint.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include "Herwig++/Decay/DecayMatrixElement.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace std;

using namespace Herwig;

VVHardGenerator::VVHardGenerator() 
  : realMESpinCorrelations_(true),
    power_(2.0),
    preqqbar_(2.2),preqg_(16.0),pregqbar_(11.0),
    b0_((11.-2./3.*5.)/4./Constants::pi),
    LambdaQCD_(91.118*GeV*exp(-1./2./((11.-2./3.*5.)/4./Constants::pi)/0.118)),
    min_pT_(2.*GeV),
    helicityConservation_(true)
{}

void VVHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << realMESpinCorrelations_
     << alphaS_ << power_ 
     << preqqbar_ << preqg_ << pregqbar_ 
     << b0_ << ounit(LambdaQCD_,GeV) << helicityConservation_
     << ounit( min_pT_,GeV )
     << FFPvertex_ << FFWvertex_ << FFZvertex_ << WWWvertex_ << FFGvertex_;
}

void VVHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> realMESpinCorrelations_
     >> alphaS_ >> power_ 
     >> preqqbar_ >> preqg_ >> pregqbar_ 
     >> b0_ >> iunit(LambdaQCD_,GeV) >> helicityConservation_
     >> iunit( min_pT_, GeV )
     >> FFPvertex_ >> FFWvertex_ >> FFZvertex_ >> WWWvertex_ >> FFGvertex_;
}

ClassDescription<VVHardGenerator> VVHardGenerator::initVVHardGenerator;
// Definition of the static class description member.

void VVHardGenerator::Init() {

  static ClassDocumentation<VVHardGenerator> documentation
    ("There is no documentation for the VVHardGenerator class");

  static Switch<VVHardGenerator,bool> interfaceSpinCorrelations
    ("SpinCorrelations",
     "Flag to select leading order spin correlations or a "
     "calculation taking into account the real NLO effects",
     &VVHardGenerator::realMESpinCorrelations_, 1, false, false);
  static SwitchOption interfaceSpinCorrelationsLeadingOrder
    (interfaceSpinCorrelations,
     "LeadingOrder",
     "Decay bosons using a leading order 2->2 calculation of the "
     "production spin density matrix",
     0);
  static SwitchOption interfaceSpinCorrelationsRealNLO
    (interfaceSpinCorrelations,
     "RealNLO",
     "Decay bosons using a production spin density matrix which "
     "takes into account the effects of real radiation",
     1);

  static Reference<VVHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &VVHardGenerator::alphaS_, false, false, true, false, false);

  static Parameter<VVHardGenerator,double> interfacePower
    ("Power",
     "The power for the sampling of the matrix elements",
     &VVHardGenerator::power_, 2.0, 1.0, 10.0,
     false, false, Interface::limited);

  static Parameter<VVHardGenerator,double> interfacePrefactorqqbar
    ("Prefactorqqbar",
     "The prefactor for the sampling of the q qbar channel",
     &VVHardGenerator::preqqbar_, 5.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<VVHardGenerator,double> interfacePrefactorqg
    ("Prefactorqg",
     "The prefactor for the sampling of the q g channel",
     &VVHardGenerator::preqg_, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);
  
  static Parameter<VVHardGenerator,double> interfacePrefactorgqbar
    ("Prefactorgqbar",
     "The prefactor for the sampling of the g qbar channel",
     &VVHardGenerator::pregqbar_, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<VVHardGenerator, Energy> interfacepTMin
    ("minPt",
     "The pT cut on hardest emision generation"
     "2*(1-Beta)*exp(-sqr(intrinsicpT/RMS))/sqr(RMS)",
     &VVHardGenerator::min_pT_, GeV, 2.*GeV, ZERO, 100000.0*GeV,
     false, false, Interface::limited);
}

void VVHardGenerator::doinit() {
  // get the vertices we need
  // get a pointer to the standard model object in the run
  static const tcHwSMPtr hwsm
    = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if (!hwsm) throw InitException() 
	       << "missing hwsm pointer in MEPP2VVPowheg::doinit()"
	       << Exception::abortnow;
  // get pointers to all required Vertex objects
  FFPvertex_ = hwsm->vertexFFP();
  FFZvertex_ = hwsm->vertexFFZ();
  WWWvertex_ = hwsm->vertexWWW();
  FFWvertex_ = hwsm->vertexFFW();
  FFGvertex_ = hwsm->vertexFFG();
}

HardTreePtr VVHardGenerator::generateHardest(ShowerTreePtr tree) {

  // Now we want to set these data vectors according to the particles we've
  // received from the current 2->2 hard collision:
  vector<ShowerProgenitorPtr> particlesToShower;
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit )
    particlesToShower.push_back(cit->first);

  qProgenitor_  = particlesToShower[0];
  qbProgenitor_ = particlesToShower[1];
  quark_        = particlesToShower[0]->progenitor();
  antiquark_    = particlesToShower[1]->progenitor();
  qHadron_      = particlesToShower[0]->beam();
  qbHadron_     = particlesToShower[1]->beam();

  if(quark_->id()<0) {
    swap(qProgenitor_,qbProgenitor_);
    swap(quark_,antiquark_);
    swap(qHadron_,qbHadron_);
  }

  // In _our_ calculation we basically define the +z axis as being given 
  // by the direction of the incoming quark for q+qb & q+g processes and
  // the incoming gluon for g+qbar processes. So now we might need to flip
  // the beams, bjorken x values, colliding partons accordingly:
  flipped_ = quark_->momentum().z()<ZERO ? true : false;

  assert(tree->outgoingLines().size()==2);
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator
	cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit )
    particlesToShower.push_back(cit->first);

  V1_           = particlesToShower[2]->progenitor();
  V2_           = particlesToShower[3]->progenitor();
  gluon_        = getParticleData(ParticleID::g)->produceParticle();

  // Abort the run if V1_ and V2_ are not just pointers to different gauge bosons
  if(!V1_||!V2_) throw Exception() 
    << "VVHardGenerator::generateHardest()\n"
    << "one or both of the gauge boson pointers is null." << Exception::abortnow;
  if(!(abs(V1_->id())==24||V1_->id()==23)||!(abs(V2_->id())==24||V2_->id()==23))
    throw Exception() 
      << "VVHardGenerator::generateHardest()\nmisidentified gauge bosons" 
      << "V1_ = " << V1_->PDGName() << "\n"
      << "V2_ = " << V2_->PDGName() << "\n"
      << Exception::abortnow;
  
  // Order the gauge bosons in the same way as in the NLO calculation
  // (the same way as in the NLO matrix element):
  // W+(->e+,nu_e) W-(->e-,nu_ebar) (MCFM: 61 [nproc])
  if(V1_->id()==-24&&V2_->id()== 24) swap(V1_,V2_); 
  // W+/-(mu+,nu_mu / mu-,nu_mubar) Z(nu_e,nu_ebar) (MCFM: 72+77 [nproc])
  if(V1_->id()== 23&&abs(V2_->id())== 24) swap(V1_,V2_); 
  // *** N.B. ***
  // We should not have to do a swap in the ZZ case, even if the different
  // (off-shell) masses of the Z's are taken into account by generating
  // the Born variables using the WZ LO/NLO matrix element(s), because
  // those transformed matrix elements are, according to mathematica, 
  // symmetric in the momenta (and therefore also the masses) of the 2 Z's.

  // If the vector bosons have decayed already then we may want to
  // to get the children_ (and any associated photons) to correct
  // spin correlations:
  StepPtr theSubProcess = generator()->eventHandler()->currentStep();
  tPVector outgoing = theSubProcess->getFinalState();
  children_.clear();
  photons_.clear();
  if(outgoing.size()>=4) {
    for(unsigned int ix=0;ix<outgoing.size();ix++) 
      if(outgoing[ix]->parents()[0]&&
	 (abs(outgoing[ix]->parents()[0]->id())==24||
	  abs(outgoing[ix]->parents()[0]->id())==23)) {
	if(outgoing[ix]->id()!=ParticleID::gamma)
	  children_.push_back(outgoing[ix]);
	else
	  photons_.push_back(outgoing[ix]);
      };
    assert(children_.size()==4);
    if(children_[0]->parents()[0]!=children_[1]->parents()[0]) 
      swap(children_[0],children_[2]);
    if(children_[0]->parents()[0]!=children_[1]->parents()[0]) 
      swap(children_[0],children_[3]);
    if(children_[0]->parents()[0]->id()!=V1_->id()) {
      swap(children_[0],children_[2]);
      swap(children_[1],children_[3]);
    }
    if(children_[0]->id()<0) swap(children_[0],children_[1]);
    if(children_[2]->id()<0) swap(children_[2],children_[3]);
    assert(children_[0]->parents()[0]==children_[1]->parents()[0]);
    assert(children_[2]->parents()[0]==children_[3]->parents()[0]);
  }

  // Now we want to construct a bornVVKinematics object. The
  // constructor for that needs all 4 momenta, q, qbar, V1_, V2_ 
  // in that order, as well as the Bjorken xq and xqbar.

  // Get the momenta first:
  vector<Lorentz5Momentum> theBornMomenta;
  theBornMomenta.push_back(quark_->momentum());
  theBornMomenta.push_back(antiquark_->momentum());
  theBornMomenta.push_back(V1_->momentum());
  theBornMomenta.push_back(V2_->momentum());

  // N.B. if the quark_ travels in the -z direction the born kinematics object
  // will detect this and rotate all particles by pi about the y axis!

  // Leading order momentum fractions:
  tcPPtr qHadron  = generator()->currentEvent()->primaryCollision()->incoming().first;
  tcPPtr qbHadron = generator()->currentEvent()->primaryCollision()->incoming().second;
  assert(qHadron->children().size()>0&&qbHadron->children().size()>0);
  if(qHadron->children()[0]->id()<0) swap(qHadron,qbHadron);

  // quark and antiquark momentum fractions respectively
  double xa = quark_    ->momentum().z()/qHadron ->momentum().z();
  double xb = antiquark_->momentum().z()/qbHadron->momentum().z();

  // Create the object containing all 2->2 __kinematic__ information:
  B_ = bornVVKinematics(theBornMomenta,xa,xb);

  // lo_me_ is the colour & spin averaged n-body matrix element squared:
  lo_me_ = lo_me(true); 

  // Attempt to generate some radiative variables and their kinematics:
  vector<Lorentz5Momentum> theRealMomenta;
  channel_ = 999;
  if(!getEvent(theRealMomenta,channel_)) return HardTreePtr();

  // Set the maximum pT for subsequent emissions:
  pT_ < min_pT_ ? qProgenitor_ ->maximumpT(min_pT_) : qProgenitor_ ->maximumpT(pT_); 
  pT_ < min_pT_ ? qbProgenitor_->maximumpT(min_pT_) : qbProgenitor_->maximumpT(pT_); 

  // Determine whether the quark or antiquark emitted:
  fermionNumberOfMother_=0;
  if((channel_==0&&theRealMomenta[0].z()/theRealMomenta[4].z()>=ZERO)||
      channel_==2) fermionNumberOfMother_ =  1;
  else if((channel_==0&&theRealMomenta[0].z()/theRealMomenta[4].z()<ZERO)||
	   channel_==1) fermionNumberOfMother_ = -1;
  assert(fermionNumberOfMother_!=0);

  // If the quark in the original tree was travelling in the -z direction
  // then we need to unflip the event (flips are automatically carried out 
  // when the original quark travels in the in -z direction when the 
  // bornVVKinematics object is created):
  if(flipped_) 
    for(unsigned int ix=0;ix<theRealMomenta.size();ix++) 
      theRealMomenta[ix].rotateY(-Constants::pi);

  // Randomly rotate the event about the beam axis:
  double randomPhi(UseRandom::rnd()*2.*Constants::pi);
  for(unsigned int ix=0;ix<theRealMomenta.size();ix++) 
    theRealMomenta[ix].rotateZ(randomPhi);

  // Warn if momentum conservation is not obeyed:
  Lorentz5Momentum inMinusOut(theRealMomenta[0]+theRealMomenta[1]
                             -theRealMomenta[2]-theRealMomenta[3]
                             -theRealMomenta[4]);
  if(inMinusOut.t()>0.1*GeV||inMinusOut.x()>0.1*GeV||
     inMinusOut.y()>0.1*GeV||inMinusOut.z()>0.1*GeV)
    cout << "VVHardGenerator::generateHardest\n"
         << "Momentum imbalance in V1 V2 rest frame\n"
         << "P_in minus P_out = " << inMinusOut/GeV << endl;

  // From the radiative kinematics we now have to form ShowerParticle objects:
  ShowerParticlePtr p1;
  ShowerParticlePtr p2;
  ShowerParticlePtr k1(new_ptr(ShowerParticle(V1_->dataPtr(),true )));
  ShowerParticlePtr k2(new_ptr(ShowerParticle(V2_->dataPtr(),true )));
  ShowerParticlePtr k ;
  // q+qbar -> V1+V2+g
  if(channel_==0) {
    p1 = new_ptr(ShowerParticle(quark_->dataPtr()           ,false));
    p2 = new_ptr(ShowerParticle(antiquark_->dataPtr()       ,false));
    k  = new_ptr(ShowerParticle(gluon_->dataPtr()           ,true ));
  }
  // q+g -> V1+V2+q
  else if(channel_==1) {
    p1 = new_ptr(ShowerParticle(quark_->dataPtr()           ,false));
    p2 = new_ptr(ShowerParticle(gluon_->dataPtr()           ,false));
    k  = new_ptr(ShowerParticle(antiquark_->dataPtr()->CC() ,true ));
  }
  // g+qbar -> V1+V2+qbar
  else {
    p1 = new_ptr(ShowerParticle(gluon_->dataPtr()           ,false));
    p2 = new_ptr(ShowerParticle(antiquark_->dataPtr()       ,false));
    k  = new_ptr(ShowerParticle(quark_->dataPtr()->CC()     ,true ));
  }
  // Set the momenta of the ShowerParticlePtr's:
  p1->set5Momentum(theRealMomenta[0]);
  p2->set5Momentum(theRealMomenta[1]);
  k1->set5Momentum(theRealMomenta[2]);
  k2->set5Momentum(theRealMomenta[3]);
  k ->set5Momentum(theRealMomenta[4]);
  // Then construct another set of ShowerPointers that will be
  // useful in creating the nasonTree, using this information:
  ShowerParticlePtr mother;
  ShowerParticlePtr spacelikeSon;
  ShowerParticlePtr timelikeSon;
  ShowerParticlePtr spectator;
  if(fermionNumberOfMother_==1) {
    mother       = new_ptr(ShowerParticle(quark_->dataPtr()    ,false));
    spacelikeSon = p1;
    timelikeSon  = k;
    spectator    = p2;
  } else if(fermionNumberOfMother_==-1) {
    mother       = new_ptr(ShowerParticle(antiquark_->dataPtr(),false));
    spacelikeSon = p2;
    timelikeSon  = k;
    spectator    = p1;
  } else {
    throw Exception() << "VVHardGenerator::generateHardest()" 
		      << "Failed to determine whether the q or qbar branched"
		      <<  Exception::runerror;
  }
  Lorentz5Momentum motherMomentum(spacelikeSon->momentum()-timelikeSon->momentum());
  motherMomentum.rescaleMass();
  mother->set5Momentum(motherMomentum);

  // Now find the Sudakov object corresponding to the IS branching:
  BranchingList branchings = evolver()->splittingGenerator()->initialStateBranchings();
  long index = abs(mother->id());
  IdList br(3);
  if(channel_==0) {                    // q+qbar -> V1+V2+g
    br[0] =  abs(spacelikeSon->id()); // makes an id list q,q,g
    br[1] =  abs(mother->id());
    br[2] =  timelikeSon->id();
  }
  else if(channel_==1||channel_==2) {   // q+g -> V1+V2+q & g+qbar -> V1+V2+qbar
    br[0] =  spacelikeSon->id();      // makes an id list g,q,qbar
    br[1] =  abs(mother->id());
    br[2] = -abs(timelikeSon->id());
  }
  SudakovPtr theSudakov;
  for(BranchingList::const_iterator cit = branchings.lower_bound(index); 
      cit != branchings.upper_bound(index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0]==br[0]&&ids[1]==br[1]&&ids[2]==br[2]) {
      theSudakov=cit->second.first;
      break;
    }
  }
  if(!theSudakov) throw Exception() 
		    << "VVHardGenerator::generateHardest()\n" 
		    << "Can't find a Sudakov with an IdList matching: " 
		    <<  br[0] << " " << br[1] << " " << br[2] << "\n"
		    <<  Exception::runerror;


  // Create HardBranchingPtrs for the particles
  HardBranchingPtr spacelikeSonBranching =
    new_ptr(HardBranching(spacelikeSon,theSudakov  ,HardBranchingPtr()   ,HardBranching::Incoming ));
  HardBranchingPtr timelikeSonBranching =
    new_ptr(HardBranching(timelikeSon ,SudakovPtr(),spacelikeSonBranching,HardBranching::Outgoing));
  HardBranchingPtr spectatorBranching =
    new_ptr(HardBranching(spectator   ,SudakovPtr(),HardBranchingPtr()   ,HardBranching::Incoming ));
  HardBranchingPtr motherBranching =
    new_ptr(HardBranching(mother      ,SudakovPtr(),spacelikeSonBranching,HardBranching::Incoming ));
  HardBranchingPtr V1_Branching =
    new_ptr(HardBranching(k1          ,SudakovPtr(),HardBranchingPtr()   ,HardBranching::Outgoing));
  HardBranchingPtr V2_Branching =
    new_ptr(HardBranching(k2          ,SudakovPtr(),HardBranchingPtr()   ,HardBranching::Outgoing));

  // N.B. The DrellYanHardGenerator first adds the timelikeSonBranching as a child
  // child of the spacelikeSonBranching before adding the motherBranching. We do
  // it the other way round in accordance with PowhegEvolver::checkShowerMomentum.
  spacelikeSonBranching->addChild(motherBranching);
  spacelikeSonBranching->addChild(timelikeSonBranching);
  motherBranching->colourPartner(spectatorBranching);
  spectatorBranching->colourPartner(motherBranching);

  vector<HardBranchingPtr> spacelikeBranchings,hardBranchings;
  spacelikeBranchings.push_back(fermionNumberOfMother_ ==  1 ? 
				spacelikeSonBranching : spectatorBranching);
  spacelikeBranchings.push_back(fermionNumberOfMother_ == -1 ? 
				spacelikeSonBranching : spectatorBranching);
  hardBranchings.push_back(motherBranching);
  hardBranchings.push_back(spectatorBranching);
  hardBranchings.push_back(V1_Branching);
  hardBranchings.push_back(V2_Branching);

  // Recalculate the hard vertex for this event:
  // For spin correlations, if an emission occurs go calculate the relevant 
  // combination of amplitudes for the ProductionMatrixElement. 
  if(realMESpinCorrelations_) {
    // Here we reset the realVVKinematics n+1 momenta to be those
    // of the lab frame in order to calculate the spin correlations.
    // Note that these momenta are not used for anything else after
    // this.
    R_.p1r(theRealMomenta[0]);
    R_.p2r(theRealMomenta[1]);
    R_.k1r(theRealMomenta[2]);
    R_.k2r(theRealMomenta[3]);
    R_.kr (theRealMomenta[4]);
    if(channel_==0) t_u_M_R_qqb_hel_amp(R_,true);
    else if(channel_==1) t_u_M_R_qg_hel_amp (R_,true);
    else if(channel_==2) t_u_M_R_gqb_hel_amp(R_,true);
    recalculateVertex();
  }
  // Construct the HardTree object needed to perform the showers
  HardTreePtr nasonTree=new_ptr(HardTree(hardBranchings,spacelikeBranchings,
					 ShowerInteraction::QCD));

  if(nasonTree->branchings().size()!=4) throw Exception() 
         << "VVHardGenerator::generateHardest()\n" 
         << "The nasonTree has " << nasonTree->branchings().size() << "branchings\n"
	 << nasonTree << "\n" <<  Exception::runerror;
  if((motherBranching->parent()!=spacelikeSonBranching)&&
     spacelikeSonBranching->parent()!=HardBranchingPtr()&&
     spectatorBranching->parent()!=HardBranchingPtr()) throw Exception() 
         << "VVHardGenerator::generateHardest()\n" 
         << "The parent-child relationships are not valid.\n" 
	 << "motherBranching->parent()       = " << motherBranching->parent()       << "\n"
	 << "spacelikeSonBranching           = " << spacelikeSonBranching           << "\n"
	 << "spectatorBranching->parent()    = " << spectatorBranching->parent()    << "\n"
	 << "spacelikeSonBranching->parent() = " << spacelikeSonBranching->parent() << "\n"
	 <<  Exception::runerror;

  if(fermionNumberOfMother_== 1) {
    nasonTree->connect(quark_    ,motherBranching   );
    nasonTree->connect(antiquark_,spectatorBranching);
    spacelikeSonBranching->beam(qProgenitor_ ->original()->parents()[0]);
    motherBranching      ->beam(qProgenitor_ ->original()->parents()[0]);
    spectatorBranching   ->beam(qbProgenitor_->original()->parents()[0]);
  } else if(fermionNumberOfMother_==-1) {
    nasonTree->connect(antiquark_,motherBranching   );
    nasonTree->connect(quark_    ,spectatorBranching);
    spacelikeSonBranching->beam(qbProgenitor_->original()->parents()[0]);
    motherBranching      ->beam(qbProgenitor_->original()->parents()[0]);
    spectatorBranching   ->beam(qProgenitor_ ->original()->parents()[0]);
  }
//   nasonTree->connect(V1_ ,V1_Branching );
//   nasonTree->connect(V2_ ,V2_Branching );

  // This if {...} else if {...} puts the mother and spectator on the same colour
  // line. If we don't do this, then when reconstructFinalStateShower calls
  // setInitialEvolutionScales it says it failed to set the colour partners, so
  // it can't set the scale and it just forgets the emission / event. This seems
  // like an unintrusive work-around until reconstructFinalStateShower is sorted.
  ColinePtr bornColourLine=new_ptr(ColourLine());
  if(fermionNumberOfMother_== 1) {
    bornColourLine->addColoured(mother);
    bornColourLine->addAntiColoured(spectator);
  } else if(fermionNumberOfMother_==-1) {
    bornColourLine->addAntiColoured(mother);
    bornColourLine->addColoured(spectator);
  }
  ShowerParticleVector particles;
  for(set<HardBranchingPtr>::iterator cit=nasonTree->branchings().begin();
      cit!=nasonTree->branchings().end();++cit) {
    particles.push_back((*cit)->branchingParticle());
  }
  evolver()->showerModel()->partnerFinder()->
     setInitialEvolutionScales(particles,true,ShowerInteraction::QCD,true);
  // Calculate the shower variables
  evolver()->showerModel()->kinematicsReconstructor()
    ->deconstructHardJets(nasonTree,evolver(),ShowerInteraction::QCD);
  vector<bool> matched(particlesToShower.size(),false);
  for(set<HardBranchingPtr>::const_iterator cit=nasonTree->branchings().begin();
      cit!=nasonTree->branchings().end();++cit) {
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
  Boost toRest,fromRest;
  toRest   = -(mother->momentum()+spectator ->momentum()).boostVector();
  fromRest =  (quark_->momentum()+antiquark_->momentum()).boostVector();
  LorentzRotation R(toRest);
  R.boost(fromRest);
  map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator tit;
  for(tit  = tree->treelinks().begin();
      tit != tree->treelinks().end();++tit) {
    ShowerTreePtr decayTree = tit->first;
    map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
      cit = decayTree->incomingLines().begin();
    // reset the momentum of the decay particle
    Lorentz5Momentum oldMomentum = cit->first->progenitor()->momentum();
    Lorentz5Momentum newMomentum = tit->second.second->momentum();
    cit->first->progenitor()->set5Momentum(newMomentum);
    cit->first->original()  ->set5Momentum(newMomentum);
    cit->first->copy()      ->set5Momentum(newMomentum);
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
    // reset the momenta of the decay products, 
    // just use the booosted ones if no QED radiation
    if(decayTree->outgoingLines().size()==2) {
      for(cjt=decayTree->outgoingLines().begin();
	  cjt!=decayTree->outgoingLines().end();++cjt) {
	for(unsigned int ix=0;ix<children_.size();++ix) {
	  if(cjt->first->original()!=children_[ix]) continue;
	  Lorentz5Momentum newChild = R*children_[ix]->momentum();
	  cjt->first->progenitor()->set5Momentum(newChild);
	  cjt->first->original()  ->set5Momentum(newChild);
	  cjt->first->copy()      ->set5Momentum(newChild);
	}
      }
    }
    // otherwise retain the direction or one of the fermions
    // at random in the boson rest frame
    else {
      LorentzRotation boostToORF(newMomentum.findBoostToCM(),
				 newMomentum.e()/oldMomentum.mass());
      tPPtr children[2];
      if(children_[0]->parents()[0]==cit->first->original()) {
	children[0] = children_[0];
	children[1] = children_[1];
      } 
      else {
	children[0] = children_[2];
	children[1] = children_[3];
      }
      if(UseRandom::rndbool()) swap(children[0],children[1]);
      double originalTheta0 = (boostToORF*children[0]->momentum()).theta();
      double originalPhi0   = (boostToORF*children[0]->momentum()).phi();
      boostToORF.rotateZ(-originalPhi0);
      boostToORF.rotateY(-originalTheta0);
      double originalPhi1   = (boostToORF*children[1]->momentum()).phi();
      LorentzRotation boost(oldMomentum.findBoostToCM(),oldMomentum.e()/oldMomentum.mass());
      tPPtr newChildren[2];
      for(cjt=decayTree->outgoingLines().begin();
	  cjt!=decayTree->outgoingLines().end();++cjt) {
	if(cjt->first->progenitor()->id()==children[0]->id())
	  newChildren[0] = cjt->first->progenitor();
	else if(cjt->first->progenitor()->id()==children[1]->id())
	  newChildren[1] = cjt->first->progenitor();
      }
      boost.rotateZ(-(boost*newChildren[0]->momentum()).phi());
      boost.rotateY(-(boost*newChildren[0]->momentum()).theta());
      boost.rotateZ(-(boost*newChildren[1]->momentum()).phi());
      boost.rotateZ( originalPhi1);
      boost.rotateY( originalTheta0);
      boost.rotateZ( originalPhi0);
      boost.boost(-newMomentum.findBoostToCM(),
 		  newMomentum.e()/oldMomentum.mass());
      for(cjt=decayTree->outgoingLines().begin();
	  cjt!=decayTree->outgoingLines().end();++cjt) {
	Lorentz5Momentum ptemp = boost*cjt->first->progenitor()->momentum();
	cjt->first->progenitor()->set5Momentum(ptemp);
	cjt->first->original()  ->set5Momentum(ptemp);
	cjt->first->copy()      ->set5Momentum(ptemp);
      }
    }
  }
  return nasonTree;
}
   
bool VVHardGenerator::canHandle(ShowerTreePtr tree) {

  // Check for two incoming particles; sometimes there will 
  // be one incoming W or Z, which comes about because the
  // code is looking for stuff to shower in the decay. 
  if(tree->incomingLines().size()!=2) return false;

  // One should be a quark and the other an antiquark:
  vector<ShowerParticlePtr> treeParticles;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin(); cit!=tree->incomingLines().end(); ++cit)
    treeParticles.push_back(cit->first->progenitor());
  if(!(treeParticles[0]->id()>0&&treeParticles[0]->id()< 6 &&
       treeParticles[1]->id()<0&&treeParticles[1]->id()>-6) &&
     !(treeParticles[1]->id()>0&&treeParticles[1]->id()< 6 &&
       treeParticles[0]->id()<0&&treeParticles[0]->id()>-6))
    return false;

  // Check that there are two outgoing particles:
  if(tree->outgoingLines().size()!=2) return false;

  // Store all the outgoing particles in treeParticles:
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  treeParticles.resize(0);
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt)
    treeParticles.push_back(cjt->first->progenitor());

  // Check that the outgoing particles are W's or Z's:
  if(!(abs(treeParticles[0]->id())==23||abs(treeParticles[0]->id())==24)) return false;
  if(!(abs(treeParticles[1]->id())==23||abs(treeParticles[1]->id())==24)) return false;

  // OK we should be able to proceed with these incoming / outgoing lines:
  return true;
}

double VVHardGenerator::getResult(int channel, realVVKinematics R, Energy pT) {

  // This routine should return the integrand of the exact Sudakov form factor,
  // defined precisely as
  // \Delta(pT) = exp[ - \int_{pT}^{pTmax} dpT dYk d\phi/2pi * getResult(...) ]
 
  // Get pi for the prefactor:
  using Constants::pi;

  // Get the VV invariant mass^2:
  Energy2 p2 = B_.sb();

  // Get the momentum fractions for the n+1 body event:
  double x1 = R.x1r();
  double x2 = R.x2r();

  // Reject the event if the x1 and x2 values are outside the phase space:
  if(x1<0.||x1>1.||x2<0.||x2>1.||x1*x2<p2/sqr(generator()->maximumCMEnergy())) return 0.;

  // Get the momentum fractions for the n body event:
  double x1b = B_.x1b();
  double x2b = B_.x2b();

  // Get the mandelstam variables needed to calculate the n+1 body matrix element:
  Energy2 s  = R.sr() ;
  Energy2 t_u_MR_o_MB;
  double  lo_lumi, nlo_lumi;

  // The luminosity function for the leading order n-body process:
  lo_lumi = qHadron_ ->pdf()->xfx(qHadron_ ,quark_    ->dataPtr(),PDFScale_,x1b)
          * qbHadron_->pdf()->xfx(qbHadron_,antiquark_->dataPtr(),PDFScale_,x2b);

  // Now we calculate the luminosity functions (product of the two pdfs) for the
  // real emission process(es) and also their matrix elements:
  // q + qbar -> V + V + g
  if(channel==0) {
    nlo_lumi    = qHadron_ ->pdf()->xfx(qHadron_ ,quark_    ->dataPtr()         ,PDFScale_,x1)
                * qbHadron_->pdf()->xfx(qbHadron_,antiquark_->dataPtr()         ,PDFScale_,x2);
    t_u_MR_o_MB =  t_u_M_R_qqb_hel_amp(R,false)/lo_me_;
  }
  // q + g -> V + V + q
  else if(channel==1) {
    nlo_lumi    = qHadron_ ->pdf()->xfx(qHadron_ ,quark_->dataPtr()             ,PDFScale_,x1)
                * qbHadron_->pdf()->xfx(qbHadron_,getParticleData(ParticleID::g),PDFScale_,x2);
    t_u_MR_o_MB =  t_u_M_R_qg_hel_amp(R,false)/lo_me_;
  }
  // g + qbar -> V + V + qbar
  else {
    nlo_lumi    = qHadron_ ->pdf()->xfx(qHadron_ ,getParticleData(ParticleID::g),PDFScale_,x1)
                * qbHadron_->pdf()->xfx(qbHadron_,antiquark_->dataPtr()         ,PDFScale_,x2);
    t_u_MR_o_MB =  t_u_M_R_gqb_hel_amp(R,false)/lo_me_;
  }  

  // Multiply ratio of the real emission matrix elements to the Born matrix element
  // by the ratio of the pdfs for the real emission and born processes to get theWeight
  if(lo_lumi<=0.||nlo_lumi<=0.) 
    return 0.;
  else 
    return t_u_MR_o_MB * ( nlo_lumi/lo_lumi * p2/s ) 
                       * sqr(p2/s)/8./pi/pi 
                       / pT / p2
                       * GeV;
} 

bool VVHardGenerator::getEvent(vector<Lorentz5Momentum> & theRealMomenta, 
				     unsigned int & channel){

  // Invariant mass of the colliding hadrons:
  Energy2 S  = sqr(generator()->maximumCMEnergy());

  // Born variables which are preserved (mass and rapidity of the diboson system):
  Energy2 p2 = B_.sb();
  double  Yb = B_.Yb();

  // Born variables which are not preserved but are needed (the momentum fractions):
  double x1b(B_.x1b()), x2b(B_.x2b());
  double x12b(x1b*x1b), x22b(x2b*x2b);

  // Maximum jet pT (half of the hadronic C.O.M. energy. N.B. this is overestimated a lot):
  Energy starting_pT = sqrt(S)/2.;

  // Initialize the pT_ *integration limit* i.e. the pT of the generated emission:
  pT_ = ZERO;

  // The pT *integration variable*:
  Energy pT;

  // The x_1 & x_2 momentum fractions corresponding to incoming momenta p1 & p2:
  double x1_(-999.), x2_(-999.);
  double x1 (-999.), x2 (-999.);

  // The jet rapidity *integration variable* and its limits:
  double Yk, minYk(-8.0), maxYk(8.0);

  // The theta2 integration variable (the azimuthal angle of the gluon w.r.t
  // V1 in the V1 & V2 rest frame:
  double theta2;

  // The realVVKinematics object corresponding to the current integration
  // set of integration variables:
  realVVKinematics R;

  // The veto algorithm rejection weight and a corresponding flag:
  double rejectionWeight;
  bool   rejectEmission ;

  // Initialize the flag indicating the selected radiation channel:
  channel=999;

  // Some product of constants used for the crude distribution:
  double a(0.);

  for(int j=0;j<3;j++) {
    pT=starting_pT;
    a =(maxYk-minYk)*prefactor_[j]/2./b0_;
    do {
      // Generate next pT according to exp[- \int^{pTold}_{pT} dpT a*(power-1)/(pT^power)]
      // 	pT = GeV/pow(  pow(GeV/pT,power_-1.) - log(UseRandom::rnd())/a
      // 		    ,  1./(power_-1.) );
      // Generate next pT according to exp[- \int^{pTold}_{pT} dpT alpha1loop*prefactor/pT ]
      pT = LambdaQCD_*exp( 0.5*exp( log(log(sqr(pT/LambdaQCD_)))+log(UseRandom::rnd())/a ) );
      // Generate rapidity of the jet:
      Yk = minYk + UseRandom::rnd()*(maxYk - minYk);
      // Generate the theta2 radiative variable:
      theta2 = UseRandom::rnd() * 2.*Constants::pi;
      // eT of the diboson system:
      Energy eT = sqrt(pT*pT+p2);
      // Calculate the eT and then solve for x_{\oplus} & x_{\ominus}:
      x1 = (pT*exp( Yk)+eT*exp( Yb))/sqrt(S);
      x2 = (pT*exp(-Yk)+eT*exp(-Yb))/sqrt(S);
      // Calculate the xr radiative variable:
      double xr(p2/(x1*x2*S));
      // Then use this to calculate the y radiative variable:
      double y(-((xr+1.)/(xr-1.))*(xr*sqr(x1/x1b)-1.)/(xr*sqr(x1/x1b)+1.));
      // The y value above should equal the one commented out below this line:
      // double y( ((xr+1.)/(xr-1.))*(xr*sqr(x2/x2b)-1.)/(xr*sqr(x2/x2b)+1.));
      // Now we get the lower limit on the x integration, xbar:
      double omy(1.-y), opy(1.+y);
      double xbar1 = 2.*opy*x12b/(sqrt(sqr(1.+x12b)*sqr(omy)+16.*y*x12b)+omy*(1.-x1b)*(1.+x1b));
      double xbar2 = 2.*omy*x22b/(sqrt(sqr(1.+x22b)*sqr(opy)-16.*y*x22b)+opy*(1.-x2b)*(1.+x2b));
      double xbar  = max(xbar1,xbar2);
      // Now we can calculate xtilde:
      double xt    = (xr-xbar)/(1.-xbar);
      // Finally we can make the realVVKinematics object:
      R = realVVKinematics(B_,xt,y,theta2);
      // The next thing we have to do is set the QCD, EW and PDF scales using R:
      setTheScales(pT);
      // ... and so calculate rejection weight:
      rejectionWeight = getResult(j,R,pT);
      // If generating according to exp[- \int^{pTold}_{pT} dpT a*(power-1)/(pT^power)]
      // rejectionWeight/= alphaS_->overestimateValue()*prefactor_[j]*pow(GeV/pT,power_);
      // If generating according to exp[- \int^{pTold}_{pT} dpT alpha1loop*prefactor/pT ]
      rejectionWeight/= 1./b0_/log(sqr(pT/LambdaQCD_))*prefactor_[j]*GeV/pT;
      rejectEmission  = UseRandom::rnd()>rejectionWeight;
      // The event is a no-emission event if pT goes past min_pT_ - basically set to 
      // outside the histogram bounds (hopefully histogram objects just ignore it then).
      if(pT<min_pT_) {
	pT=ZERO;
	rejectEmission = false;
      }
      if(rejectionWeight>1.) {
	ostringstream stream;
	stream << "VVHardGenerator::getEvent weight for channel " << j
	       << " is greater than one: " << rejectionWeight << endl;
	generator()->logWarning( Exception(stream.str(), Exception::warning) );
      }
    }
    while(rejectEmission);
    // set pT of emission etc
    if(pT>pT_) {
      channel = j;
      pT_ = pT;
      Yk_ = Yk;
      R_  = R ;
      x1_ = x1;
      x2_ = x2;
    }
  }
  // Was this an (overall) no emission event?
  if(pT_<min_pT_) { 
    pT_ = ZERO;
    channel = 3;
  }
  if(channel==3) return false;
  if(channel>3) throw Exception() 
	       << "VVHardGenerator::getEvent() channel = " << channel
	       << "   pT = " << pT/GeV << "   pT_ = " << pT_/GeV
 	       << Exception::abortnow;

  // Work out the momenta in the lab frame, reserving the mass and rapidity 
  // of the VV system:
  LorentzRotation yzRotation;
  yzRotation.setRotateX(-atan2(pT_/GeV,sqrt(p2)/GeV));
  LorentzRotation boostFrompTisZero;
  boostFrompTisZero.setBoostY(-pT_/sqrt(p2+pT_*pT_));
  LorentzRotation boostFromYisZero;
  boostFromYisZero.setBoostZ(tanh(Yb));

  theRealMomenta.resize(5);
  theRealMomenta[0] = Lorentz5Momentum(ZERO,ZERO, x1_*sqrt(S)/2., x1_*sqrt(S)/2.,ZERO);
  theRealMomenta[1] = Lorentz5Momentum(ZERO,ZERO,-x2_*sqrt(S)/2., x2_*sqrt(S)/2.,ZERO);
  theRealMomenta[2] = boostFromYisZero*boostFrompTisZero*yzRotation*(R_.k1r());
  theRealMomenta[3] = boostFromYisZero*boostFrompTisZero*yzRotation*(R_.k2r());
  theRealMomenta[4] = Lorentz5Momentum(ZERO, pT_,  pT_*sinh(Yk_),  pT_*cosh(Yk_),ZERO);

  return true;
}

void VVHardGenerator::doinitrun() {
  // insert the different prefactors in the vector for easy look up
  prefactor_.push_back(preqqbar_);
  prefactor_.push_back(preqg_);
  prefactor_.push_back(pregqbar_);
  HardestEmissionGenerator::doinitrun();
}

void VVHardGenerator::setTheScales(Energy pT) {
  // Work out the scales we want to use in the matrix elements and the pdfs:
  // Scale for alpha_S: pT^2 of the diboson system.
  QCDScale_ = max(pT*pT,sqr(min_pT_));
  // Scale for real emission PDF: 
  // pT^2+mVV^2 - as mcfm does in the case of a single W/Z boson).
  // Energy2 PDFScale_ = max(R.pT2_in_lab(),sqr(min_pT_))+R.s2r();
  // pT^2 - as advocated by Nason & Ridolfi for ZZ production & Alioli et al for gg->h: 
  PDFScale_ = max(pT*pT,sqr(min_pT_));
  // Scale of electroweak vertices: mVV^2 the invariant mass of the diboson system.
  //  EWScale_  = B_.sb();
  // ... And this choice is more like what can be seen in mcatnlo_vbmain.f (weird).
  EWScale_ = 0.5*(B_.k12b()+B_.k22b());

  return;
}

/***************************************************************************/
// This is identical to the code in the Powheg matrix element. It should
// equal t_u_M_R_qqb in there, which is supposed to be the real emission ME
// times tk*uk.
Energy2 VVHardGenerator::t_u_M_R_qqb_hel_amp(realVVKinematics R, bool getMatrix) const {
  using namespace ThePEG::Helicity;

  ProductionMatrixElement qqb_hel_amps(PDT::Spin1Half,PDT::Spin1Half,
				       PDT::Spin1    ,PDT::Spin1    ,
				       PDT::Spin1);

  double sum_hel_amps_sqr(0.);

  tcPDPtr p1data(quark_->dataPtr());
  tcPDPtr p2data(antiquark_->dataPtr());
  tcPDPtr k1data(V1_->dataPtr());
  tcPDPtr k2data(V2_->dataPtr());
  tcPDPtr kdata(getParticleData(ParticleID::g));
  if(k1data->id()==-24&&k2data->id()==24) swap(k1data,k2data);

  SpinorWaveFunction qSpinor(R.p1r(),p1data,incoming);
  SpinorBarWaveFunction qbSpinor(R.p2r(),p2data,incoming);
  vector<SpinorWaveFunction> q;
  vector<SpinorBarWaveFunction> qb;
  for(unsigned int ix=0;ix<2;ix++) {
    qSpinor.reset(ix);
    qbSpinor.reset(ix);
    q.push_back(qSpinor);
    qb.push_back(qbSpinor);
  }

  VectorWaveFunction v1Polarization(R.k1r(),k1data,outgoing);
  VectorWaveFunction v2Polarization(R.k2r(),k2data,outgoing);
  vector<VectorWaveFunction> v1;
  vector<VectorWaveFunction> v2;
  for(unsigned int ix=0;ix<3;ix++) {
    v1Polarization.reset(ix);
    v2Polarization.reset(ix);
    v1.push_back(v1Polarization);
    v2.push_back(v2Polarization);
  }

  VectorWaveFunction gPolarization(R.kr(),kdata,outgoing);
  vector<VectorWaveFunction> g;
  for(unsigned int ix=0;ix<3;ix+=2) {
    gPolarization.reset(ix);
    g.push_back(gPolarization);
  }

  AbstractFFVVertexPtr ffg  = FFGvertex_;
  AbstractFFVVertexPtr ffv1 = k1data->id()==23 ? FFZvertex_ : FFWvertex_;
  AbstractFFVVertexPtr ffv2 = k2data->id()==23 ? FFZvertex_ : FFWvertex_;

  // Collecting information for intermediate fermions
  vector<tcPDPtr> tc;
  if(abs(k1data->id())==24&&abs(k2data->id())==24) {
    if(abs(p1data->id())%2==0)
      for(int ix=0;ix<3;++ix) tc.push_back(getParticleData(1+2*ix));
    else
      for(int ix=0;ix<3;++ix) tc.push_back(getParticleData(2+2*ix));
  }
  else if(k1data->id()==23&&k2data->id()==23)      tc.push_back(p1data);
  else if(abs(k1data->id())==24&&k2data->id()==23) tc.push_back(p2data);

  // Loop over helicities summing the relevant diagrams
  for(unsigned int p1hel=0;p1hel<2;++p1hel) {
    for(unsigned int p2hel=0;p2hel<2;++p2hel) {
      for(unsigned int khel=0;khel<2;++khel) {
	SpinorWaveFunction    p1_k  = ffg->evaluate(QCDScale_,5,p1data,q[p1hel],g[khel]);
	SpinorBarWaveFunction p2_k  = ffg->evaluate(QCDScale_,5,p2data,qb[p2hel],g[khel]);
	for(unsigned int k1hel=0;k1hel<3;++k1hel) {
	  for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	    // If helicity is exactly conserved (massless quarks) skip if p1hel=p2hel
	    // but if the production ME is required first fill it with (0.,0.).
 	    if((p1hel==p2hel)&&helicityConservation_) {
	      if(getMatrix) {
		if(khel==0)
		  qqb_hel_amps(p1hel,p2hel,k1hel,k2hel,0) = Complex(0.,0.);
		else
		  qqb_hel_amps(p1hel,p2hel,k1hel,k2hel,2) = Complex(0.,0.);
	      }
 	      continue;
 	    }
	    vector<Complex> diagrams;
	    // Get all t-channel diagram contributions
	    tcPDPtr intermediate_t;
	    for(unsigned int ix=0;ix<tc.size();ix++) {
	      intermediate_t = (!(k1data->id()==24&&k2data->id()==-24)) ? p2data : tc[ix];
	    // Note: choosing 5 as the second argument ffvX_->evaluate() sets
	    // option 5 in thepeg/Helicity/Vertex/VertexBase.cc, which makes
	    // the (fermion) propagator denominator massless: 1/p^2.
	    // If W+Z / W-Z calculate the two V+jet-like s-channel diagrams
	      SpinorWaveFunction    p1_v1 = ffv1->evaluate(EWScale_,5,intermediate_t,q[p1hel],v1[k1hel]);
	      SpinorBarWaveFunction p2_v2 = ffv2->evaluate(EWScale_,5,intermediate_t,qb[p2hel],v2[k2hel]);
	      // First calculate all the off-shell fermion currents
	      // Now calculate the 6 t-channel diagrams
	      // q+qb->g+v1+v2, q+qb->v1+g+v2, q+qb->v1+v2+g
	      if(!((k1data->id()==24&&k2data->id()==-24)&&(abs(p1data->id())%2==1))) {
		diagrams.push_back(ffv1->evaluate(EWScale_,p1_k,p2_v2,v1[k1hel]));
		diagrams.push_back(ffg->evaluate(QCDScale_,p1_v1,p2_v2,g[khel]));
		diagrams.push_back(ffv2->evaluate(EWScale_,p1_v1,p2_k,v2[k2hel]));
	      }
	      intermediate_t = (!(k1data->id()==24&&k2data->id()==-24)) ? p1data : tc[ix];
	      SpinorWaveFunction    p1_v2 = ffv2->evaluate(EWScale_,5,intermediate_t,q[p1hel],v2[k2hel]);
	      SpinorBarWaveFunction p2_v1 = ffv1->evaluate(EWScale_,5,intermediate_t,qb[p2hel],v1[k1hel]);
	      // q+qb->g+v2+v1, q+qb->v2+g+v1, q+qb->v2+v1+g
	      if(!((k1data->id()==24&&k2data->id()==-24)&&(abs(p1data->id())%2==0))) {
		diagrams.push_back(ffv2->evaluate(EWScale_,p1_k,p2_v1,v2[k2hel]));
		diagrams.push_back(ffg->evaluate(QCDScale_,p1_v2,p2_v1,g[khel]));
		diagrams.push_back(ffv1->evaluate(EWScale_,p1_v2,p2_k,v1[k1hel]));
	      }
	    }
	    // Note: choosing 3 as the second argument in WWWvertex_->evaluate() 
	    // sets option 3 in thepeg/Helicity/Vertex/VertexBase.cc , which 
	    // means the propagator does not contain a width factor (which is 
	    // good re. gauge invariance). 
	    // If W+Z / W-Z calculate the two V+jet-like s-channel diagrams
	    if(abs(k1data->id())==24&&k2data->id()==23) {
	      // The off-shell s-channel boson current
	      VectorWaveFunction k1_k2 = 
		WWWvertex_->evaluate(EWScale_,3,k1data->CC(),v2[k2hel],v1[k1hel]);
	      // q+qb->g+v1*->g+v1+v2, q+qb->v1*+g->v1+v2+g
	      diagrams.push_back(ffv1->evaluate(EWScale_,p1_k,qb[p2hel],k1_k2));
	      diagrams.push_back(ffv1->evaluate(EWScale_,q[p1hel],p2_k,k1_k2));
	    }
	    // If W+W- calculate the four V+jet-like s-channel diagrams
	    if((k1data->id()==24&&k2data->id()==-24)&&(p1data->id()==-p2data->id())) {
	      // The off-shell s-channel boson current
	      VectorWaveFunction k1_k2;
	      // q+qb->g+Z0*->g+v1+v2,q+qb->Z0*+g->v1+v2+g,
	      tcPDPtr Z0    = getParticleData(ParticleID::Z0);
	      k1_k2 = WWWvertex_->evaluate(EWScale_,3,Z0,v2[k2hel],v1[k1hel]);
	      diagrams.push_back(FFZvertex_->evaluate(EWScale_,p1_k,qb[p2hel],k1_k2));
	      diagrams.push_back(FFZvertex_->evaluate(EWScale_,q[p1hel],p2_k,k1_k2));
	      // q+qb->g+gamma*->g+v1+v2,q+qb->gamma*+g->v1+v2+g,
	      tcPDPtr gamma = getParticleData(ParticleID::gamma);
	      k1_k2 = WWWvertex_->evaluate(EWScale_,3,gamma,v2[k2hel],v1[k1hel]);
	      diagrams.push_back(FFPvertex_->evaluate(EWScale_,p1_k,qb[p2hel],k1_k2));
	      diagrams.push_back(FFPvertex_->evaluate(EWScale_,q[p1hel],p2_k,k1_k2));
	    }
	    // Add up all diagrams to get the total amplitude:
	    Complex hel_amp(0.);
	    for(unsigned int ix=0;ix<diagrams.size();ix++) hel_amp += diagrams[ix];
	    // If we need to fill the production ME we do it here:
 	    if(getMatrix) {
	      if(khel==0)
		qqb_hel_amps(p1hel,p2hel,k1hel,k2hel,0) = hel_amp;
	      else
		qqb_hel_amps(p1hel,p2hel,k1hel,k2hel,2) = hel_amp;
	    }
	    sum_hel_amps_sqr += norm(hel_amp);
	  }
	}
      }
    }
  }

  // Fill up the remaining bits of the production ME, corresponding 
  // to longitudinal gluon polarization, with (0.,0.).
  if(getMatrix) {
    for(unsigned int p1hel=0;p1hel<2;++p1hel) {
      for(unsigned int p2hel=0;p2hel<2;++p2hel) {
	for(unsigned int k1hel=0;k1hel<3;++k1hel) {
	  for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	    qqb_hel_amps(p1hel,p2hel,k1hel,k2hel,1) = Complex(0.,0.);
	  }
	}
      }
    }
  }

  // Calculate the production density matrix:
  if(getMatrix) {
    for(unsigned int k1hel=0;k1hel<3;++k1hel) {
      for(unsigned int k1helpr=0;k1helpr<3;++k1helpr) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	  for(unsigned int k2helpr=0;k2helpr<3;++k2helpr) {
	    Complex theElement(0.,0.);
	    // For each k1hel, k1helpr, k2hel, k2helpr sum over fermion and gluon spins...
	    for(unsigned int p1hel=0;p1hel<2;++p1hel) {
	      for(unsigned int p2hel=0;p2hel<2;++p2hel) {
		for(unsigned int khel=0;khel<3;khel+=2) {
		  theElement += qqb_hel_amps(p1hel,p2hel,k1hel  ,k2hel  ,khel)
		          *conj(qqb_hel_amps(p1hel,p2hel,k1helpr,k2helpr,khel));
		}
	      }
	    }
	    // ...and then set the production matrix element to the sum:
	    productionMatrix_[k1hel][k1helpr][k2hel][k2helpr] = theElement;
	  }
	}
      }
    }
  }

  // Spin and colour averaging factors = 1/4 * CF * 1/3 = 1/9
  sum_hel_amps_sqr /= 9.;

  // Symmetry factor for identical Z bosons in the final state 
  if(k1data->id()==23&&k2data->id()==23) sum_hel_amps_sqr /= 2.;

  return sum_hel_amps_sqr*R.tkr()*R.ukr()*UnitRemoval::InvE2;
}

/***************************************************************************/
// This is identical to the code in the Powheg matrix element. It should
// equal t_u_M_R_qg in there, which is supposed to be the real emission ME
// times tk*uk.
Energy2 VVHardGenerator::t_u_M_R_qg_hel_amp(realVVKinematics R, bool getMatrix) const {
  using namespace ThePEG::Helicity;

  ProductionMatrixElement qg_hel_amps(PDT::Spin1Half,PDT::Spin1,
				      PDT::Spin1,PDT::Spin1,
				      PDT::Spin1Half);
  
  double sum_hel_amps_sqr(0.);

  tcPDPtr p1data(quark_->dataPtr());
  tcPDPtr p2data(getParticleData(ParticleID::g));
  tcPDPtr k1data(V1_->dataPtr());
  tcPDPtr k2data(V2_->dataPtr());
  tcPDPtr kdata (antiquark_->dataPtr()->CC());
  if(k1data->id()==-24&&k2data->id()==24) swap(k1data,k2data);

  SpinorWaveFunction qinSpinor(R.p1r(),p1data,incoming);
  SpinorBarWaveFunction qoutSpinor(R.kr(),kdata,outgoing);
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qout;
  for(unsigned int ix=0;ix<2;ix++) {
    qinSpinor.reset(ix);
    qoutSpinor.reset(ix);
    qin.push_back(qinSpinor);
    qout.push_back(qoutSpinor);
  }

  VectorWaveFunction v1Polarization(R.k1r(),k1data,outgoing);
  VectorWaveFunction v2Polarization(R.k2r(),k2data,outgoing);
  vector<VectorWaveFunction> v1;
  vector<VectorWaveFunction> v2;
  for(unsigned int ix=0;ix<3;ix++) {
    v1Polarization.reset(ix);
    v2Polarization.reset(ix);
    v1.push_back(v1Polarization);
    v2.push_back(v2Polarization);
  }

  VectorWaveFunction gPolarization(R.p2r(),p2data,incoming);
  vector<VectorWaveFunction> g;
  for(unsigned int ix=0;ix<3;ix+=2) {
    gPolarization.reset(ix);
    g.push_back(gPolarization);
  }

  AbstractFFVVertexPtr ffg  = FFGvertex_;
  AbstractFFVVertexPtr ffv1 = k1data->id()==23 ? FFZvertex_ : FFWvertex_;
  AbstractFFVVertexPtr ffv2 = k2data->id()==23 ? FFZvertex_ : FFWvertex_;

  // Collecting information for intermediate fermions
  vector<tcPDPtr> tc;
  if(abs(k1data->id())==24&&abs(k2data->id())==24) {
    if(abs(p1data->id())%2==0)
      for(int ix=0;ix<3;++ix) tc.push_back(getParticleData(1+2*ix));
    else
      for(int ix=0;ix<3;++ix) tc.push_back(getParticleData(2+2*ix));
  }
  else if(k1data->id()==23&&k2data->id()==23)      tc.push_back(p1data);
  else if(abs(k1data->id())==24&&k2data->id()==23) tc.push_back(kdata->CC());

  // Loop over helicities summing the relevant diagrams
  for(unsigned int p1hel=0;p1hel<2;++p1hel) {
    for(unsigned int p2hel=0;p2hel<2;++p2hel) {
      for(unsigned int khel=0;khel<2;++khel) {
	SpinorWaveFunction    p1_p2 = ffg->evaluate(QCDScale_,5,p1data,qin[p1hel],g[p2hel]);
	SpinorBarWaveFunction p2_k  = ffg->evaluate(QCDScale_,5,kdata->CC(),qout[khel],g[p2hel]);
	for(unsigned int k1hel=0;k1hel<3;++k1hel) {
	  for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	    // If helicity is exactly conserved (massless quarks) skip if p1hel!=khel
	    // but if the production ME is required first fill it with (0.,0.).
	    if((p1hel!=khel)&&helicityConservation_) {
	      if(getMatrix) {
		if(p2hel==0)
		  qg_hel_amps(p1hel,0,k1hel,k2hel,khel) = Complex(0.,0.);
		else
		  qg_hel_amps(p1hel,2,k1hel,k2hel,khel) = Complex(0.,0.);
	      }
	      continue;
	    }
	    vector<Complex> diagrams;
	    // Get all t-channel diagram contributions
	    tcPDPtr intermediate_q;
	    for(unsigned int ix=0;ix<tc.size();ix++) {
	      intermediate_q = (!(k1data->id()==24&&k2data->id()==-24)) ? antiquark_->dataPtr() : tc[ix];
	      SpinorWaveFunction    p1_v1 = ffv1->evaluate(EWScale_,5,intermediate_q,qin[p1hel],v1[k1hel]);
	      SpinorBarWaveFunction k_v2  = ffv2->evaluate(EWScale_,5,intermediate_q,qout[khel],v2[k2hel]);
	      // First calculate all the off-shell fermion currents
	      // Now calculate the 6 abelian diagrams
	      // q+g->v1+v2+q with 2 t-channel propagators, 1 s- and 1 t-channel and 2 t-channel ones.
	      if(!((k1data->id()==24&&k2data->id()==-24)&&(abs(p1data->id())%2==1))) {
		diagrams.push_back(ffv2->evaluate(EWScale_,p1_v1,p2_k,v2[k2hel]));
		diagrams.push_back(ffg->evaluate(QCDScale_,p1_v1,k_v2,g[p2hel]));
		diagrams.push_back(ffv1->evaluate(EWScale_,p1_p2,k_v2,v1[k1hel]));
	      }
	      intermediate_q = (!(k1data->id()==24&&k2data->id()==-24)) ? p1data : tc[ix];
	      SpinorWaveFunction    p1_v2 = ffv2->evaluate(EWScale_,5,intermediate_q,qin[p1hel],v2[k2hel]);
              SpinorBarWaveFunction k_v1  = ffv1->evaluate(EWScale_,5,intermediate_q,qout[khel],v1[k1hel]);
	      // q+g->v2+v1+q, with 2 t-channel propagators, 1 s- and 1 t-channel and 2 t-channel ones.
	      if(!((k1data->id()==24&&k2data->id()==-24)&&(abs(p1data->id())%2==0))) {
		diagrams.push_back(ffv1->evaluate(EWScale_,p1_v2,p2_k,v1[k1hel]));
		diagrams.push_back(ffg->evaluate(QCDScale_,p1_v2,k_v1,g[p2hel]));
		diagrams.push_back(ffv2->evaluate(EWScale_,p1_p2,k_v1,v2[k2hel]));
	      }
	    }
	    // Note: choosing 3 as the second argument in WWWvertex_->evaluate() 
	    // sets option 3 in thepeg/Helicity/Vertex/VertexBase.cc , which 
	    // means the propagator does not contain a width factor (which is 
	    // good re. gauge invariance). 
	    // If W+Z / W-Z calculate the two V+jet-like s-channel diagrams
	    if(abs(k1data->id())==24&&k2data->id()==23) {
	      // The off-shell s-channel boson current
	      VectorWaveFunction k1_k2 = 
		WWWvertex_->evaluate(EWScale_,3,k1data->CC(),v2[k2hel],v1[k1hel]);
	      // q+qb->g+v1*->g+v1+v2, q+qb->v1*+g->v1+v2+g
	      diagrams.push_back(ffv1->evaluate(EWScale_,p1_p2,qout[khel],k1_k2));
	      diagrams.push_back(ffv1->evaluate(EWScale_,qin[p1hel],p2_k,k1_k2));
	    }
	    // If W+W- calculate the four V+jet-like s-channel diagrams
	    if((k1data->id()==24&&k2data->id()==-24)&&(p1data->id()==kdata->id())) {
	      // The off-shell s-channel boson current
	      VectorWaveFunction k1_k2;
	      // q+qb->g+Z0*->g+v1+v2,q+qb->Z0*+g->v1+v2+g,
	      tcPDPtr Z0    = getParticleData(ParticleID::Z0);
	      k1_k2 = WWWvertex_->evaluate(EWScale_,3,Z0,v2[k2hel],v1[k1hel]);
	      diagrams.push_back(FFZvertex_->evaluate(EWScale_,p1_p2,qout[khel],k1_k2));
	      diagrams.push_back(FFZvertex_->evaluate(EWScale_,qin[p1hel],p2_k,k1_k2));
	      // q+qb->g+gamma*->g+v1+v2,q+qb->gamma*+g->v1+v2+g,
	      tcPDPtr gamma = getParticleData(ParticleID::gamma);
	      k1_k2 = WWWvertex_->evaluate(EWScale_,3,gamma,v2[k2hel],v1[k1hel]);
	      diagrams.push_back(FFPvertex_->evaluate(EWScale_,p1_p2,qout[khel],k1_k2));
	      diagrams.push_back(FFPvertex_->evaluate(EWScale_,qin[p1hel],p2_k,k1_k2));
	    }
	    // Add up all diagrams to get the total amplitude:
	    Complex hel_amp(0.);
	    for(unsigned int ix=0;ix<diagrams.size();ix++) hel_amp += diagrams[ix];
	    // If we need to fill the production ME we do it here:
 	    if(getMatrix) {
	      if(p2hel==0)
		qg_hel_amps(p1hel,0,k1hel,k2hel,khel) = hel_amp;
	      else
		qg_hel_amps(p1hel,2,k1hel,k2hel,khel) = hel_amp;
	    }
	    sum_hel_amps_sqr += norm(hel_amp);
	  }
	}
      }
    }
  }
  
  // Fill up the remaining bits of the production ME, corresponding 
  // to longitudinal gluon polarization, with (0.,0.).
  if(getMatrix) {
    for(unsigned int p1hel=0;p1hel<2;++p1hel) {
      for(unsigned int k1hel=0;k1hel<3;++k1hel) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	  for(unsigned int khel=0;khel<2;++khel) {
	    qg_hel_amps(p1hel,1,k1hel,k2hel,khel) = Complex(0.,0.);
	  }
	}
      }
    }
  }

  // Calculate the production density matrix:
  if(getMatrix) {
    for(unsigned int k1hel=0;k1hel<3;++k1hel) {
      for(unsigned int k1helpr=0;k1helpr<3;++k1helpr) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	  for(unsigned int k2helpr=0;k2helpr<3;++k2helpr) {
	    Complex theElement(0.,0.);
	    // For each k1hel, k1helpr, k2hel, k2helpr sum over fermion and gluon spins...
	    for(unsigned int p1hel=0;p1hel<2;++p1hel) {
	      for(unsigned int p2hel=0;p2hel<3;p2hel+=2) {
		for(unsigned int khel=0;khel<2;++khel) {
		  theElement += qg_hel_amps(p1hel,p2hel,k1hel  ,k2hel  ,khel)
		          *conj(qg_hel_amps(p1hel,p2hel,k1helpr,k2helpr,khel));
		}
	      }
	    }
	    // ...and then set the production matrix element to the sum:
	    productionMatrix_[k1hel][k1helpr][k2hel][k2helpr] = theElement;
	  }
	}
      }
    }
  }

  // Spin and colour averaging factors = 1/4 * TR * 1/3 = 1/24
  sum_hel_amps_sqr /= 24.;

  // Symmetry factor for identical Z bosons in the final state 
  if(k1data->id()==23&&k2data->id()==23) sum_hel_amps_sqr /= 2.;

  return sum_hel_amps_sqr*R.tkr()*R.ukr()*UnitRemoval::InvE2;
}

/***************************************************************************/
// This is identical to the code in the Powheg matrix element. It should
// equal t_u_M_R_gqb in there, which is supposed to be the real emission ME
// times tk*uk.
Energy2 VVHardGenerator::t_u_M_R_gqb_hel_amp(realVVKinematics R, bool getMatrix) const {
  using namespace ThePEG::Helicity;

  ProductionMatrixElement gqb_hel_amps(PDT::Spin1,PDT::Spin1Half,
				       PDT::Spin1,PDT::Spin1,
				       PDT::Spin1Half);
  
  double sum_hel_amps_sqr(0.);

  tcPDPtr p1data(getParticleData(ParticleID::g));
  tcPDPtr p2data(antiquark_->dataPtr());
  tcPDPtr k1data(V1_->dataPtr());
  tcPDPtr k2data(V2_->dataPtr());
  tcPDPtr kdata (quark_->dataPtr()->CC());
  if(k1data->id()==-24&&k2data->id()==24) swap(k1data,k2data);

  SpinorBarWaveFunction qbinSpinor(R.p2r(),p2data,incoming);
  SpinorWaveFunction qboutSpinor(R.kr(),kdata,outgoing);
  vector<SpinorBarWaveFunction> qbin;
  vector<SpinorWaveFunction> qbout;
  for(unsigned int ix=0;ix<2;ix++) {
    qbinSpinor.reset(ix);
    qboutSpinor.reset(ix);
    qbin.push_back(qbinSpinor);
    qbout.push_back(qboutSpinor);
  }

  VectorWaveFunction v1Polarization(R.k1r(),k1data,outgoing);
  VectorWaveFunction v2Polarization(R.k2r(),k2data,outgoing);
  vector<VectorWaveFunction> v1;
  vector<VectorWaveFunction> v2;
  for(unsigned int ix=0;ix<3;ix++) {
    v1Polarization.reset(ix);
    v2Polarization.reset(ix);
    v1.push_back(v1Polarization);
    v2.push_back(v2Polarization);
  }

  VectorWaveFunction gPolarization(R.p1r(),p1data,incoming);
  vector<VectorWaveFunction> g;
  for(unsigned int ix=0;ix<3;ix+=2) {
    gPolarization.reset(ix);
    g.push_back(gPolarization);
  }

  AbstractFFVVertexPtr ffg  = FFGvertex_;
  AbstractFFVVertexPtr ffv1 = k1data->id()==23 ? FFZvertex_ : FFWvertex_;
  AbstractFFVVertexPtr ffv2 = k2data->id()==23 ? FFZvertex_ : FFWvertex_;

  // Collecting information for intermediate fermions
  vector<tcPDPtr> tc;
  if(abs(k1data->id())==24&&abs(k2data->id())==24) {
    if(abs(p2data->id())%2==0)
      for(int ix=0;ix<3;++ix) tc.push_back(getParticleData(1+2*ix));
    else
      for(int ix=0;ix<3;++ix) tc.push_back(getParticleData(2+2*ix));
  }
  else if(k1data->id()==23&&k2data->id()==23)      tc.push_back(p2data);
  else if(abs(k1data->id())==24&&k2data->id()==23) tc.push_back(kdata->CC());

  // Loop over helicities summing the relevant diagrams
  for(unsigned int p1hel=0;p1hel<2;++p1hel) {
    for(unsigned int p2hel=0;p2hel<2;++p2hel) {
      for(unsigned int khel=0;khel<2;++khel) {
	SpinorBarWaveFunction p1_p2 = ffg->evaluate(QCDScale_,5,p2data,qbin[p2hel],g[p1hel]);
	SpinorWaveFunction    p1_k  = ffg->evaluate(QCDScale_,5,kdata->CC(),qbout[khel],g[p1hel]);
	for(unsigned int k1hel=0;k1hel<3;++k1hel) {
	  for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	    // If helicity is exactly conserved (massless quarks) skip if p2hel!=khel
	    // but if the production ME is required first fill it with (0.,0.).
	    if((p2hel!=khel)&&helicityConservation_) {
	      if(getMatrix) {
		if(p1hel==0)
		  gqb_hel_amps(0,p2hel,k1hel,k2hel,khel) = Complex(0.,0.);
		else
		  gqb_hel_amps(2,p2hel,k1hel,k2hel,khel) = Complex(0.,0.);
	      }
	      continue;
	    }
	    vector<Complex> diagrams;
	    // Get all t-channel diagram contributions
	    tcPDPtr intermediate_q;
	    for(unsigned int ix=0;ix<tc.size();ix++) {
	      intermediate_q = (!(k1data->id()==24&&k2data->id()==-24)) ? quark_->dataPtr() : tc[ix];
	      SpinorBarWaveFunction p2_v1 = ffv1->evaluate(EWScale_,5,intermediate_q,qbin[p2hel],v1[k1hel]);
	      SpinorWaveFunction    k_v2  = ffv2->evaluate(EWScale_,5,intermediate_q,qbout[khel],v2[k2hel]);
	      // First calculate all the off-shell fermion currents
	      // Now calculate the 6 abelian diagrams q+g->v1+v2+q 
              // with 2 t-channel propagators, 1 s- and 1 t-channel 
              // and 2 t-channel ones.
	      if(!((k1data->id()==24&&k2data->id()==-24)&&(abs(p2data->id())%2==0))) {
		diagrams.push_back(ffv2->evaluate(EWScale_,p1_k,p2_v1,v2[k2hel]));
		diagrams.push_back(ffg->evaluate(QCDScale_,k_v2,p2_v1,g[p1hel]));
		diagrams.push_back(ffv1->evaluate(EWScale_,k_v2,p1_p2,v1[k1hel]));
	      }
	      intermediate_q = (!(k1data->id()==24&&k2data->id()==-24)) ? p2data : tc[ix];
	      SpinorBarWaveFunction p2_v2 = ffv2->evaluate(EWScale_,5,intermediate_q,qbin[p2hel],v2[k2hel]);
              SpinorWaveFunction    k_v1  = ffv1->evaluate(EWScale_,5,intermediate_q,qbout[khel],v1[k1hel]);
	      // q+g->v2+v1+q, with 2 t-channel propagators, 1 s- and 1 t-channel and 2 t-channel ones.
	      if(!((k1data->id()==24&&k2data->id()==-24)&&(abs(p2data->id())%2==1))) {
		diagrams.push_back(ffv1->evaluate(EWScale_,p1_k,p2_v2,v1[k1hel]));
		diagrams.push_back(ffg->evaluate(QCDScale_,k_v1,p2_v2,g[p1hel]));
		diagrams.push_back(ffv2->evaluate(EWScale_,k_v1,p1_p2,v2[k2hel]));
	      }
	    }
	    // Note: choosing 3 as the second argument in WWWvertex_->evaluate() 
	    // sets option 3 in thepeg/Helicity/Vertex/VertexBase.cc , which 
	    // means the propagator does not contain a width factor (which is 
	    // good re. gauge invariance). 
	    // If W+Z / W-Z calculate the two V+jet-like s-channel diagrams
	    if(abs(k1data->id())==24&&k2data->id()==23) {
	      // The off-shell s-channel boson current
	      VectorWaveFunction k1_k2 = 
		WWWvertex_->evaluate(EWScale_,3,k1data->CC(),v2[k2hel],v1[k1hel]);
	      // q+qb->g+v1*->g+v1+v2, q+qb->v1*+g->v1+v2+g
	      diagrams.push_back(ffv1->evaluate(EWScale_,qbout[khel],p1_p2,k1_k2));
	      diagrams.push_back(ffv1->evaluate(EWScale_,p1_k,qbin[p2hel],k1_k2));
	    }
	    // If W+W- calculate the four V+jet-like s-channel diagrams
	    if((k1data->id()==24&&k2data->id()==-24)&&(p2data->id()==kdata->id())) {
	      // The off-shell s-channel boson current
	      VectorWaveFunction k1_k2;
	      // q+qb->g+Z0*->g+v1+v2,q+qb->Z0*+g->v1+v2+g,
	      tcPDPtr Z0    = getParticleData(ParticleID::Z0);
	      k1_k2 = WWWvertex_->evaluate(EWScale_,3,Z0,v2[k2hel],v1[k1hel]);
	      diagrams.push_back(FFZvertex_->evaluate(EWScale_,qbout[khel],p1_p2,k1_k2));
	      diagrams.push_back(FFZvertex_->evaluate(EWScale_,p1_k,qbin[p2hel],k1_k2));
	      // q+qb->g+gamma*->g+v1+v2,q+qb->gamma*+g->v1+v2+g,
	      tcPDPtr gamma = getParticleData(ParticleID::gamma);
	      k1_k2 = WWWvertex_->evaluate(EWScale_,3,gamma,v2[k2hel],v1[k1hel]);
	      diagrams.push_back(FFPvertex_->evaluate(EWScale_,qbout[khel],p1_p2,k1_k2));
	      diagrams.push_back(FFPvertex_->evaluate(EWScale_,p1_k,qbin[p2hel],k1_k2));
	    }
	    // Add up all diagrams to get the total amplitude:
	    Complex hel_amp(0.);
	    for(unsigned int ix=0;ix<diagrams.size();ix++) hel_amp += diagrams[ix];
	    // If we need to fill the production ME we do it here:
 	    if(getMatrix) {
	      if(p1hel==0)
		gqb_hel_amps(0,p2hel,k1hel,k2hel,khel) = hel_amp;
	      else
		gqb_hel_amps(2,p2hel,k1hel,k2hel,khel) = hel_amp;
	    }
	    sum_hel_amps_sqr += norm(hel_amp);
	  }
	}
      }
    }
  }
  
  // Fill up the remaining bits of the production ME, corresponding 
  // to longitudinal gluon polarization, with (0.,0.).
  if(getMatrix) {
    for(unsigned int p2hel=0;p2hel<2;++p2hel) {
      for(unsigned int k1hel=0;k1hel<3;++k1hel) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	  for(unsigned int khel=0;khel<2;++khel) {
	    gqb_hel_amps(1,p2hel,k1hel,k2hel,khel) = Complex(0.,0.);
	  }
	}
      }
    }
  }

  // Calculate the production density matrix:
  if(getMatrix) {
    for(unsigned int k1hel=0;k1hel<3;++k1hel) {
      for(unsigned int k1helpr=0;k1helpr<3;++k1helpr) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	  for(unsigned int k2helpr=0;k2helpr<3;++k2helpr) {
	    Complex theElement(0.,0.);
	    // For each k1hel, k1helpr, k2hel, k2helpr sum over fermion and gluon spins...
	    for(unsigned int p1hel=0;p1hel<3;p1hel+=2) {
	      for(unsigned int p2hel=0;p2hel<2;++p2hel) {
		for(unsigned int khel=0;khel<2;++khel) {
		  theElement += gqb_hel_amps(p1hel,p2hel,k1hel  ,k2hel  ,khel)
		          *conj(gqb_hel_amps(p1hel,p2hel,k1helpr,k2helpr,khel));
		}
	      }
	    }
	    // ...and then set the production matrix element to the sum:
	    productionMatrix_[k1hel][k1helpr][k2hel][k2helpr] = theElement;
	  }
	}
      }
    }
  }

  // Spin and colour averaging factors = 1/4 * TR * 1/3 = 1/24
  sum_hel_amps_sqr /= 24.;

  // Symmetry factor for identical Z bosons in the final state 
  if(k1data->id()==23&&k2data->id()==23) sum_hel_amps_sqr /= 2.;

  return sum_hel_amps_sqr*R.tkr()*R.ukr()*UnitRemoval::InvE2;
}


/***************************************************************************/
// This returns exactly the same value as lo_me2_ when you put it in MEPP2VVPowheg.cc
double VVHardGenerator::lo_me(bool getMatrix) const {
  using namespace ThePEG::Helicity;

  ProductionMatrixElement lo_hel_amps(PDT::Spin1Half,PDT::Spin1Half,
				      PDT::Spin1    ,PDT::Spin1);

  double sum_hel_amps_sqr(0.);

  tcPDPtr p1data(quark_->dataPtr());
  tcPDPtr p2data(antiquark_->dataPtr());
  tcPDPtr k1data(V1_->dataPtr());
  tcPDPtr k2data(V2_->dataPtr());
  if(k1data->id()==-24&&k2data->id()==24) swap(k1data,k2data); // Should never actually occur.

  // If you want to reproduce the spin correlations of MEPP2VV
  // you should evaluate this ME using the lab frame momenta
  // instead of the bornVVKinematics ones (partonic C.O.M. frame).
  SpinorWaveFunction qSpinor;
  SpinorBarWaveFunction qbSpinor;
  if(!getMatrix) {
    qSpinor=SpinorWaveFunction(B_.p1b(),p1data,incoming);
    qbSpinor=SpinorBarWaveFunction(B_.p2b(),p2data,incoming);
  } else {
    qSpinor=SpinorWaveFunction(quark_->momentum(),p1data,incoming);
    qbSpinor=SpinorBarWaveFunction(antiquark_->momentum(),p2data,incoming);
  }
  vector<SpinorWaveFunction> q;
  vector<SpinorBarWaveFunction> qb;
  for(unsigned int ix=0;ix<2;ix++) {
    qSpinor.reset(ix);
    qbSpinor.reset(ix);
    q.push_back(qSpinor);
    qb.push_back(qbSpinor);
  }

  // If you want to reproduce the spin correlations of MEPP2VV
  // you should evaluate this ME using the lab frame momenta
  // instead of the bornVVKinematics ones (partonic C.O.M. frame).
  VectorWaveFunction v1Polarization;
  VectorWaveFunction v2Polarization;
  if(!getMatrix) {
    v1Polarization=VectorWaveFunction(B_.k1b(),k1data,outgoing);
    v2Polarization=VectorWaveFunction(B_.k2b(),k2data,outgoing);
  } else {
    v1Polarization=VectorWaveFunction(V1_->momentum(),k1data,outgoing);
    v2Polarization=VectorWaveFunction(V2_->momentum(),k2data,outgoing);
  }
  vector<VectorWaveFunction> v1;
  vector<VectorWaveFunction> v2;
  for(unsigned int ix=0;ix<3;ix++) {
    v1Polarization.reset(ix);
    v2Polarization.reset(ix);
    v1.push_back(v1Polarization);
    v2.push_back(v2Polarization);
  }

  AbstractFFVVertexPtr ffv1 = k1data->id()==23 ? FFZvertex_ : FFWvertex_;
  AbstractFFVVertexPtr ffv2 = k2data->id()==23 ? FFZvertex_ : FFWvertex_;

  // Collecting information for intermediate fermions
  vector<tcPDPtr> tc;
  if(abs(k1data->id())==24&&abs(k2data->id())==24) {
    if(abs(p1data->id())%2==0)
      for(int ix=0;ix<3;++ix) tc.push_back(getParticleData(1+2*ix));
    else
      for(int ix=0;ix<3;++ix) tc.push_back(getParticleData(2+2*ix));
  }
  else if(k1data->id()==23&&k2data->id()==23)      tc.push_back(p1data);
  else if(abs(k1data->id())==24&&k2data->id()==23) tc.push_back(p2data);

  // Loop over helicities summing the relevant diagrams
  for(unsigned int p1hel=0;p1hel<2;++p1hel) {
    for(unsigned int p2hel=0;p2hel<2;++p2hel) {
      for(unsigned int k1hel=0;k1hel<3;++k1hel) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	  if((p1hel==p2hel)&&helicityConservation_) {
	    lo_hel_amps(p1hel,p2hel,k1hel,k2hel) = Complex(0.,0.);
	    continue;
	  }
	  vector<Complex> diagrams;
	  // Get all t-channel diagram contributions
	  tcPDPtr intermediate_t;
	  for(unsigned int ix=0;ix<tc.size();ix++) {
	    intermediate_t = (!(k1data->id()==24&&k2data->id()==-24)) ? p2data : tc[ix];
	    SpinorWaveFunction    p1_v1 = ffv1->evaluate(EWScale_,5,intermediate_t,q[p1hel],v1[k1hel]);
	    // First calculate all the off-shell fermion currents
	    // Now calculate the 6 t-channel diagrams
	    // q+qb->v1+v2
	    if(!((k1data->id()==24&&k2data->id()==-24)&&(abs(p1data->id())%2==1)))
	      diagrams.push_back(ffv2->evaluate(EWScale_,p1_v1,qb[p2hel],v2[k2hel]));
	    intermediate_t = (!(k1data->id()==24&&k2data->id()==-24)) ? p1data : tc[ix];
	    SpinorWaveFunction    p1_v2 = ffv2->evaluate(EWScale_,5,intermediate_t,q[p1hel],v2[k2hel]);
	    // q+qb->v2+v1
	    if(!((k1data->id()==24&&k2data->id()==-24)&&(abs(p1data->id())%2==0)))
	      diagrams.push_back(ffv1->evaluate(EWScale_,p1_v2,qb[p2hel],v1[k1hel]));
	  }
	  // Note: choosing 3 as the second argument in WWWvertex_->evaluate() 
	  // sets option 3 in thepeg/Helicity/Vertex/VertexBase.cc , which 
	  // means the propagator does not contain a width factor (which is 
	  // good re. gauge invariance). 
	  // If W+Z / W-Z calculate the two V+jet-like s-channel diagrams
	  if(abs(k1data->id())==24&&k2data->id()==23) {
	    // The off-shell s-channel boson current
	    VectorWaveFunction k1_k2 = 
	      WWWvertex_->evaluate(EWScale_,3,k1data->CC(),v2[k2hel],v1[k1hel]);
	    // q+qb->v1*->v1+v2
	    diagrams.push_back(ffv1->evaluate(EWScale_,q[p1hel],qb[p2hel],k1_k2));
	  }
	  // If W+W- calculate the four V+jet-like s-channel diagrams
	  if((k1data->id()==24&&k2data->id()==-24)&&(p1data->id()==-p2data->id())) {
	    // The off-shell s-channel boson current
	    VectorWaveFunction k1_k2;
	    // q+qb->Z0*->v1+v2
	    tcPDPtr Z0    = getParticleData(ParticleID::Z0);
	    k1_k2 = WWWvertex_->evaluate(EWScale_,3,Z0,v2[k2hel],v1[k1hel]);
	    diagrams.push_back(FFZvertex_->evaluate(EWScale_,q[p1hel],qb[p2hel],k1_k2));
	    // q+qb->gamma*->v1+v2
	    tcPDPtr gamma = getParticleData(ParticleID::gamma);
	    k1_k2 = WWWvertex_->evaluate(EWScale_,3,gamma,v2[k2hel],v1[k1hel]);
	    diagrams.push_back(FFPvertex_->evaluate(EWScale_,q[p1hel],qb[p2hel],k1_k2));
	  }
	  // Add up all diagrams to get the total amplitude:
	  Complex hel_amp(0.);
	  for(unsigned int ix=0;ix<diagrams.size();ix++) hel_amp += diagrams[ix];
	  // If we need to fill the production ME we do it here:
	  if(getMatrix) lo_hel_amps(p1hel,p2hel,k1hel,k2hel) = hel_amp;
	  sum_hel_amps_sqr += norm(hel_amp);
	}
      }
    }
  }

  // Calculate the production density matrix:
  if(getMatrix) {
    for(unsigned int k1hel=0;k1hel<3;++k1hel) {
      for(unsigned int k1helpr=0;k1helpr<3;++k1helpr) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	  for(unsigned int k2helpr=0;k2helpr<3;++k2helpr) {
	    Complex theElement(0.,0.);
	    // For each k1hel, k1helpr, k2hel, k2helpr sum over the fermion spins...
	    for(unsigned int p1hel=0;p1hel<2;++p1hel) {
	      for(unsigned int p2hel=0;p2hel<2;++p2hel) {
		if((p1hel==p2hel)&&helicityConservation_) continue;
		theElement += lo_hel_amps(p1hel,p2hel,k1hel  ,k2hel  )
		        *conj(lo_hel_amps(p1hel,p2hel,k1helpr,k2helpr));
	      }
	    }
	    // ...and then set the production matrix element to the sum:
	    productionMatrix_[k1hel][k1helpr][k2hel][k2helpr] = theElement;
	  }
	}
      }
    }
  }

  // Spin and colour averaging factors = 1/4 * 1/3 = 1/12
  sum_hel_amps_sqr /= 12.;

  // Symmetry factor for identical Z bosons in the final state 
  if(k1data->id()==23&&k2data->id()==23) sum_hel_amps_sqr /= 2.;

  return sum_hel_amps_sqr;
}

/***************************************************************************/
// This member selects a [2-body] decay mode and assigns children to the
// vector bosons with momenta which are isotropic in their rest frames.
bool VVHardGenerator::isotropicDecayer() {
  using namespace ThePEG::Helicity;

  // Generate the children's momenta isotropically in
  // the rest frames of V1 and V2:
  double cth,phi;
  // First V1's children:
  cth = UseRandom::rnd()*2.-1.;
  phi = UseRandom::rnd()*2.*Constants::pi;
  Energy m1(V1_->momentum().m());
  Energy m3(children_[0]->data().constituentMass());
  Energy m4(children_[1]->data().constituentMass());
  Energy p34(triangleFn(sqr(m1),sqr(m3),sqr(m4))
	     /2./m1);
  if(isnan(p34/GeV)||cth>1.||cth<-1.) return false;
  Energy pT34(p34*sqrt(1.-cth)*sqrt(1.+cth));
  Lorentz5Momentum k3(pT34*sin(phi),pT34*cos(phi),p34 *cth,
		      sqrt(p34*p34+sqr(m3)),m3);
  Lorentz5Momentum k4(-k3);
  k4.setE(sqrt(p34*p34+sqr(m4)));
  k4.setTau(m4);
  Boost boostToV1RF(R_.k1r().boostVector());
  k3.boost(boostToV1RF);
  k3.rescaleRho();
  k4.boost(boostToV1RF);
  k4.rescaleRho();

  // Second V2's children:
  cth = UseRandom::rnd()*2.-1.;
  phi = UseRandom::rnd()*2.*Constants::pi;
  Energy m2(V2_->momentum().m());
  Energy m5(children_[2]->data().constituentMass());
  Energy m6(children_[3]->data().constituentMass());
  Energy p56(triangleFn(sqr(m2),sqr(m5),sqr(m6))
	     /2./m2);
  if(isnan(p56/GeV)||cth>1.||cth<-1.) return false;
  Energy pT56(p56*sqrt(1.-cth)*sqrt(1.+cth));
  Lorentz5Momentum k5(pT56*sin(phi),pT56*cos(phi),p56*cth,
		      sqrt(p56*p56+sqr(m5)),m5);
  Lorentz5Momentum k6(-k5);
  k6.setE(sqrt(p56*p56+sqr(m6)));
  k6.setTau(m6);
  Boost boostToV2RF(R_.k2r().boostVector());
  k5.boost(boostToV2RF);
  k5.rescaleRho();
  k6.boost(boostToV2RF);
  k6.rescaleRho();

  // Assign the momenta to the children:
  children_[0]->set5Momentum(k3);
  children_[1]->set5Momentum(k4);
  children_[2]->set5Momentum(k5);
  children_[3]->set5Momentum(k6);
  return true;

}

// Override 2->2 production matrix here:
void VVHardGenerator::recalculateVertex() {

  // Zero the squared amplitude; this equals sum_hel_amps_sqr if all
  // is working as it should:
  Complex productionMatrix2(0.,0.);
  for(unsigned int k1hel=0;k1hel<3;++k1hel)
    for(unsigned int k2hel=0;k2hel<3;++k2hel)
      productionMatrix2 += productionMatrix_[k1hel][k1hel][k2hel][k2hel];

  // Get the vector wavefunctions:
  VectorWaveFunction v1Polarization;
  VectorWaveFunction v2Polarization;
  v1Polarization=VectorWaveFunction(R_.k1r(),V1_->dataPtr(),outgoing);
  v2Polarization=VectorWaveFunction(R_.k2r(),V2_->dataPtr(),outgoing);
  vector<VectorWaveFunction> v1;
  vector<VectorWaveFunction> v2;
  for(unsigned int ix=0;ix<3;ix++) {
    v1Polarization.reset(ix);
    v2Polarization.reset(ix);
    v1.push_back(v1Polarization);
    v2.push_back(v2Polarization);
  }

  AbstractFFVVertexPtr ffv1 = V1_->id()==23 ? FFZvertex_ : FFWvertex_;
  AbstractFFVVertexPtr ffv2 = V2_->id()==23 ? FFZvertex_ : FFWvertex_;

  bool vetoed(true);
  while(vetoed) {
    // Decay the bosons isotropically in their rest frames:
    isotropicDecayer();

    // Get the spinor wavefunctions:
    SpinorWaveFunction k3Spinor(children_[0]->momentum(),children_[0]->dataPtr(),outgoing);
    SpinorBarWaveFunction k4Spinor(children_[1]->momentum(),children_[1]->dataPtr(),outgoing);
    SpinorWaveFunction k5Spinor(children_[2]->momentum(),children_[2]->dataPtr(),outgoing);
    SpinorBarWaveFunction k6Spinor(children_[3]->momentum(),children_[3]->dataPtr(),outgoing);
    vector<SpinorWaveFunction> k3,k5;
    vector<SpinorBarWaveFunction> k4,k6;
    for(unsigned int ix=0;ix<2;ix++) {
      k3Spinor.reset(ix);
      k4Spinor.reset(ix);
      k3.push_back(k3Spinor);
      k4.push_back(k4Spinor);
      k5Spinor.reset(ix);
      k6Spinor.reset(ix);
      k5.push_back(k5Spinor);
      k6.push_back(k6Spinor);
    }

    DecayMatrixElement decayAmps(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half);

    for(unsigned int k1hel=0;k1hel<3;++k1hel) {
      for(unsigned int k3hel=0;k3hel<2;++k3hel) {
	for(unsigned int k4hel=0;k4hel<2;++k4hel) {
	  decayAmps(k1hel,k3hel,k4hel) = 
	    ffv1->evaluate(EWScale_,k3[k3hel],k4[k4hel],v1[k1hel]);
	}
      }
    }
    Complex V1decayMatrix[3][3];
    for(unsigned int k1hel=0;k1hel<3;++k1hel) {
      for(unsigned int k1helpr=0;k1helpr<3;++k1helpr) {
	Complex theElement(0.,0.);
	for(unsigned int k3hel=0;k3hel<2;++k3hel) {
	  for(unsigned int k4hel=0;k4hel<2;++k4hel) {
	    theElement += decayAmps(k1hel,k3hel,k4hel)
	            *conj(decayAmps(k1helpr,k3hel,k4hel));
	  }
	}
	V1decayMatrix[k1hel][k1helpr] = theElement;
      }
    }
    Complex V1decayMatrix2(0.,0.);
    for(unsigned int k1hel=0;k1hel<3;++k1hel) V1decayMatrix2 += V1decayMatrix[k1hel][k1hel];

    for(unsigned int k2hel=0;k2hel<3;++k2hel) {
      for(unsigned int k5hel=0;k5hel<2;++k5hel) {
	for(unsigned int k6hel=0;k6hel<2;++k6hel) {
	  decayAmps(k2hel,k5hel,k6hel) = 
	    ffv2->evaluate(EWScale_,k5[k5hel],k6[k6hel],v2[k2hel]);
	}
      }
    }
    Complex V2decayMatrix[3][3];
    for(unsigned int k2hel=0;k2hel<3;++k2hel) {
      for(unsigned int k2helpr=0;k2helpr<3;++k2helpr) {
	Complex theElement(0.,0.);
	for(unsigned int k5hel=0;k5hel<2;++k5hel) {
	  for(unsigned int k6hel=0;k6hel<2;++k6hel) {
	    theElement += decayAmps(k2hel,k5hel,k6hel)
	            *conj(decayAmps(k2helpr,k5hel,k6hel));
	  }
	}
	V2decayMatrix[k2hel][k2helpr] = theElement;
      }
    }
    Complex V2decayMatrix2(0.,0.);
    for(unsigned int k2hel=0;k2hel<3;++k2hel) V2decayMatrix2 += V2decayMatrix[k2hel][k2hel];

    // Contract the production matrix and the decay matrices:
    Complex meTimesV1V2denominators(0.,0.);
    for(unsigned int k1hel=0;k1hel<3;++k1hel) {
      for(unsigned int k1helpr=0;k1helpr<3;++k1helpr) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	  for(unsigned int k2helpr=0;k2helpr<3;++k2helpr) {
	    meTimesV1V2denominators += 
	      productionMatrix_[k1hel][k1helpr][k2hel][k2helpr]
	      *V1decayMatrix[k1hel][k1helpr]
	      *V2decayMatrix[k2hel][k2helpr];
	  }
	}
      }
    }

    if(imag(meTimesV1V2denominators)/real(meTimesV1V2denominators)>1.e-7) 
      cout << "VVHardGenerator warning\n" 
	   << "the matrix element's imaginary part is large " 
	   << meTimesV1V2denominators << endl;  
    if(imag(productionMatrix2)/real(productionMatrix2)>1.e-7) 
      cout << "VVHardGenerator warning\n" 
	   << "the production matrix element's imaginary part is large " 
	   << productionMatrix2 << endl;  
    if(imag(V1decayMatrix2)/real(V1decayMatrix2)>1.e-7) 
      cout << "VVHardGenerator warning\n" 
	   << "the V1 decay matrix element's imaginary part is large " 
	   << V1decayMatrix2 << endl;  
    if(imag(V2decayMatrix2)/real(V2decayMatrix2)>1.e-7) 
      cout << "VVHardGenerator warning\n" 
	   << "the V2 decay matrix element's imaginary part is large " 
	   << V2decayMatrix2 << endl;  

    // Need branching ratio at least in here I would think --->
    double decayWeight( real(meTimesV1V2denominators)
		      / real(productionMatrix2*V1decayMatrix2*V2decayMatrix2));
    if(decayWeight>1.)
      cout << "VVHardGenerator::recalculateVertex decayWeight > 1., decayWeight = "
	   << decayWeight << endl;
    if(decayWeight<0.)
      cout << "VVHardGenerator::recalculateVertex decayWeight < 0., decayWeight = "
	   << decayWeight << endl;
    if(UseRandom::rnd()<decayWeight) vetoed = false;
    else vetoed = true;
  }

  return;

}

Energy2 VVHardGenerator::triangleFn(Energy2 m12,Energy2 m22, Energy2 m32) {
  Energy4 lambda2(m12*m12+m22*m22+m32*m32-2.*m12*m22-2.*m12*m32-2.*m22*m32);
  if(lambda2>=ZERO) {
    return sqrt(lambda2);
  } else {
    generator()->log() 
      << "VVHardGenerator::triangleFn "
      << "kinematic instability, imaginary triangle function\n";
    return -999999.*GeV2;
  }
}
