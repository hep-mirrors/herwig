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
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"

using namespace std;

using namespace Herwig;

VVHardGenerator::VVHardGenerator() 
  : power_(2.0),
    preqqbar_(30.0),preqg_(4.0),pregqbar_(4.0),
    min_pT_(2.*GeV)
{}

void VVHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << alphaS_ << power_ 
     << preqqbar_ << preqg_ << pregqbar_ 
     << ounit( min_pT_,GeV )
     << FFPvertex_ << FFWvertex_ << FFZvertex_ << WWWvertex_ << FFGvertex_;
}

void VVHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> alphaS_ >> power_ 
     >> preqqbar_ >> preqg_ >> pregqbar_ 
     >> iunit( min_pT_, GeV )
     >> FFPvertex_ >> FFWvertex_ >> FFZvertex_ >> WWWvertex_ >> FFGvertex_;
}

ClassDescription<VVHardGenerator> VVHardGenerator::initVVHardGenerator;
// Definition of the static class description member.

void VVHardGenerator::Init() {

  static ClassDocumentation<VVHardGenerator> documentation
    ("There is no documentation for the VVHardGenerator class");

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
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit )
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

  // Reset the gauge bosons to null pointers (do not remove!):
  V1_ = PPtr();
  V2_ = PPtr();

  // Get the gauge bosons:
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    PPtr currentBoson;
    currentBoson = cjt->first->copy();
    if(!V1_&&!V2_)                    V1_ = currentBoson;
    if( V1_&&!V2_&&V1_!=currentBoson) V2_ = currentBoson;
  }
  gluon_ = getParticleData(ParticleID::g);

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
  lo_me_ = lo_me();

  // Attempt to generate some radiative variables and their kinematics:
  vector<Lorentz5Momentum> theRealMomenta;
  int channel(-999);
  if(!getEvent(theRealMomenta,channel)) return HardTreePtr();

  // Set the maximum pT for subsequent emissions:
  pT_ < min_pT_ ? qProgenitor_ ->maximumpT(min_pT_) : qProgenitor_ ->maximumpT(pT_); 
  pT_ < min_pT_ ? qbProgenitor_->maximumpT(min_pT_) : qbProgenitor_->maximumpT(pT_); 

  // Determine whether the quark or antiquark emitted:
  int fermionNumberOfMother(0);
  if((channel==0&&theRealMomenta[0].z()/theRealMomenta[4].rapidity()>ZERO)||
      channel==2) fermionNumberOfMother =  1;
  else if((channel==0&&theRealMomenta[0].z()/theRealMomenta[4].rapidity()<ZERO)||
	   channel==1) fermionNumberOfMother = -1;
  assert(fermionNumberOfMother!=0);

  // If the quark in the original tree was travelling in the -z direction
  // then we need to unflip the event (flips are automatically carried out 
  // when the original quark travels in the in -z direction when the 
  // bornVVKinematics object is created):
  if(flipped_) 
    for(unsigned int ix=0;ix<theRealMomenta.size();ix++) 
      theRealMomenta[ix].rotateY(-Constants::pi);

  // Randomly rotate the event about the beam axis:
  for(unsigned int ix=0;ix<theRealMomenta.size();ix++) 
    theRealMomenta[ix].rotateZ(UseRandom::rnd()*2.*Constants::pi);

  // From the radiative kinematics we now have to form ShowerParticle objects:
  ShowerParticlePtr p1;
  ShowerParticlePtr p2;
  ShowerParticlePtr k1(new_ptr(ShowerParticle(V1_->dataPtr(),true )));
  ShowerParticlePtr k2(new_ptr(ShowerParticle(V2_->dataPtr(),true )));
  ShowerParticlePtr k ;
  // q+qbar -> V1+V2+g
  if(channel==0) {
    p1 = new_ptr(ShowerParticle(quark_->dataPtr()           ,false));
    p2 = new_ptr(ShowerParticle(antiquark_->dataPtr()       ,false));
    k  = new_ptr(ShowerParticle(gluon_                      ,true ));
  }
  // q+g -> V1+V2+q
  else if(channel==1) {
    p1 = new_ptr(ShowerParticle(quark_->dataPtr()           ,false));
    p2 = new_ptr(ShowerParticle(gluon_                      ,false));
    k  = new_ptr(ShowerParticle(antiquark_->dataPtr()->CC() ,true ));
  }
  // g+qbar -> V1+V2+qbar
  else {
    p1 = new_ptr(ShowerParticle(gluon_                      ,false));
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
  if(fermionNumberOfMother==1) {
    mother       = new_ptr(ShowerParticle(quark_->dataPtr()    ,false));
    spacelikeSon = p1;
    timelikeSon  = k;
    spectator    = p2;
  } else if(fermionNumberOfMother==-1) {
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
  if(channel==0) {                    // q+qbar -> V1+V2+g
    br[0] =  abs(spacelikeSon->id()); // makes an id list q,q,g
    br[1] =  abs(mother->id());
    br[2] =  timelikeSon->id();
  }
  else if(channel==1||channel==2) {   // q+g -> V1+V2+q & g+qbar -> V1+V2+qbar
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
  spacelikeBranchings.push_back(fermionNumberOfMother ==  1 ? 
				spacelikeSonBranching : spectatorBranching);
  spacelikeBranchings.push_back(fermionNumberOfMother == -1 ? 
				spacelikeSonBranching : spectatorBranching);
  hardBranchings.push_back(motherBranching);
  hardBranchings.push_back(spectatorBranching);
  hardBranchings.push_back(V1_Branching);
  hardBranchings.push_back(V2_Branching);

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

  if(fermionNumberOfMother== 1) {
    nasonTree->connect(quark_    ,motherBranching   );
    nasonTree->connect(antiquark_,spectatorBranching);
    spacelikeSonBranching->beam(qProgenitor_ ->original()->parents()[0]);
    motherBranching      ->beam(qProgenitor_ ->original()->parents()[0]);
    spectatorBranching   ->beam(qbProgenitor_->original()->parents()[0]);
  } else if(fermionNumberOfMother==-1) {
    nasonTree->connect(antiquark_,motherBranching   );
    nasonTree->connect(quark_    ,spectatorBranching);
    spacelikeSonBranching->beam(qbProgenitor_->original()->parents()[0]);
    motherBranching      ->beam(qbProgenitor_->original()->parents()[0]);
    spectatorBranching   ->beam(qProgenitor_ ->original()->parents()[0]);
  }

  // This if {...} else if {...} puts the mother and spectator on the same colour
  // line. If we don't do this, then when reconstructFinalStateShower calls
  // setInitialEvolutionScales it says it failed to set the colour partners, so
  // it can't set the scale and it just forgets the emission / event. This seems
  // like an unintrusive work-around until reconstructFinalStateShower is sorted.
  ColinePtr bornColourLine=new_ptr(ColourLine());
  if(fermionNumberOfMother== 1) {
    bornColourLine->addColoured(mother);
    bornColourLine->addAntiColoured(spectator);
  } else if(fermionNumberOfMother==-1) {
    bornColourLine->addAntiColoured(mother);
    bornColourLine->addColoured(spectator);
  }

  // Calculate the shower variables
  evolver()->showerModel()->kinematicsReconstructor()
    ->deconstructHardJets(nasonTree,evolver(),ShowerInteraction::QCD);

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
  Energy2 tk = R.tkr();
  Energy2 uk = R.ukr();
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
    t_u_MR_o_MB =  t_u_M_R_qqb_hel_amp(R)/lo_me_;
  }
  // q + g -> V + V + q
  else if(channel==1) {
    nlo_lumi    = qHadron_ ->pdf()->xfx(qHadron_ ,quark_->dataPtr()             ,PDFScale_,x1)
                * qbHadron_->pdf()->xfx(qbHadron_,getParticleData(ParticleID::g),PDFScale_,x2);
    t_u_MR_o_MB =  t_u_M_R_qg_hel_amp(R)/lo_me_;
  }
  // g + qbar -> V + V + qbar
  else {
    nlo_lumi    = qHadron_ ->pdf()->xfx(qHadron_ ,getParticleData(ParticleID::g),PDFScale_,x1)
                * qbHadron_->pdf()->xfx(qbHadron_,antiquark_->dataPtr()         ,PDFScale_,x2);
    t_u_MR_o_MB =  t_u_M_R_gqb_hel_amp(R)/lo_me_;
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
				     int & channel){

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

  // The pT *integration variable* and the corresponding eT of the diboson system:
  Energy pT, eT;

  // The x and y radiative variables:
  double xr, y;

  // The x_1 & x_2 momentum fractions corresponding to incoming momenta p1 & p2:
  double x1_, x2_;
  double x1 , x2 ;

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
  channel=-1;

  for(int j=0;j<1;j++) {
    pT=starting_pT;
    double a = alphaS_->overestimateValue() * prefactor_[j] * (maxYk-minYk) / (power_ - 1.);
    do {
      // Generate next pT:
      pT = GeV/pow(pow(GeV/pT,power_-1) - log(UseRandom::rnd())/a,1./(power_-1.));
      // Generate rapidity of the jet:
      Yk = minYk + UseRandom::rnd()*(maxYk - minYk);
      // Generate the theta2 radiative variable:
      theta2 = UseRandom::rnd() * 2.*Constants::pi;
      // Calculate the eT and then solve for x_{\oplus} & x_{\ominus}:
      eT = sqrt(pT*pT+p2);
      x1 = (pT*exp( Yk)+eT*exp( Yb))/sqrt(S);
      x2 = (pT*exp(-Yk)+eT*exp(-Yb))/sqrt(S);
      // Calculate the xr radiative variable:
      xr = p2/(x1*x2*S);
      // Then use this to calculate the y radiative variable:
      y  = 1.-4.*xr*pT*pT/p2/sqr(1.-xr);
      if(exp(2.*(Yk-Yb))<1.) y *= -1.; // y<0. branch
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
      rejectionWeight/= prefactor_[j]*pow(GeV/pT,power_);
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
  EWScale_  = B_.sb();

  return;
}

/***************************************************************************/
// The game here is to get this helicity amplitude squared to return all the
// same values as t_u_M_R_qqb above, TIMES a further factor tk*uk!
Energy2 VVHardGenerator::t_u_M_R_qqb_hel_amp(realVVKinematics R) const {
  using namespace ThePEG::Helicity;

//   qqb_hel_amps_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
// 					      PDT::Spin1,PDT::Spin1,
// 					      PDT::Spin1));

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
      for(unsigned int ix=0;ix<3;++ix) tc.push_back(getParticleData(1+2*ix));
    else
      for(unsigned int ix=0;ix<3;++ix) tc.push_back(getParticleData(2+2*ix));
  }
  else if(k1data->id()==23&&k2data->id()==23)      tc.push_back(p1data);
  else if(abs(k1data->id())==24&&k2data->id()==23) tc.push_back(p2data);

  // Loop over helicities summing the relevant diagrams
  for(unsigned int p1hel=0;p1hel<2;++p1hel) {
    for(unsigned int p2hel=0;p2hel<2;++p2hel) {
      for(unsigned int k1hel=0;k1hel<3;++k1hel) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	  for(unsigned int khel=0;khel<2;++khel) {
	    vector<Complex> diagrams;
	    SpinorWaveFunction    p1_k  = ffg->evaluate(QCDScale_,5,p1data,q[p1hel],g[khel]);
	    SpinorBarWaveFunction p2_k  = ffg->evaluate(QCDScale_,5,p2data,qb[p2hel],g[khel]);
	    // Get all t-channel diagram contributions
	    tcPDPtr intermediate_t;
	    for(unsigned int ix=0;ix<tc.size();ix++) {
	      intermediate_t = (!(k1data->id()==24&&k2data->id()==-24)) ? p2data : tc[ix];
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
// 	    qqb_hel_amps_(p1hel,p2hel,k1hel,k2hel,khel) = hel_amp;
	    sum_hel_amps_sqr += norm(hel_amp);
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
// The game here is to get this helicity amplitude squared to return all the
// same values as t_u_M_R_qg above, TIMES a further factor tk*uk!
Energy2 VVHardGenerator::t_u_M_R_qg_hel_amp(realVVKinematics R) const {
  using namespace ThePEG::Helicity;

//   qg_hel_amps_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
// 					     PDT::Spin1,PDT::Spin1,
// 					     PDT::Spin1Half));
  
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
      for(unsigned int ix=0;ix<3;++ix) tc.push_back(getParticleData(1+2*ix));
    else
      for(unsigned int ix=0;ix<3;++ix) tc.push_back(getParticleData(2+2*ix));
  }
  else if(k1data->id()==23&&k2data->id()==23)      tc.push_back(p1data);
  else if(abs(k1data->id())==24&&k2data->id()==23) tc.push_back(kdata->CC());

  // Loop over helicities summing the relevant diagrams
  for(unsigned int p1hel=0;p1hel<2;++p1hel) {
    for(unsigned int p2hel=0;p2hel<2;++p2hel) {
      for(unsigned int k1hel=0;k1hel<3;++k1hel) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	  for(unsigned int khel=0;khel<2;++khel) {
	    vector<Complex> diagrams;
	    SpinorWaveFunction    p1_p2 = ffg->evaluate(QCDScale_,5,p1data,qin[p1hel],g[p2hel]);
	    SpinorBarWaveFunction p2_k  = ffg->evaluate(QCDScale_,5,kdata->CC(),qout[khel],g[p2hel]);
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
// 	    qg_hel_amps_(p1hel,p2hel,k1hel,k2hel,khel) = hel_amp;
	    sum_hel_amps_sqr += norm(hel_amp);
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
// The game here is to get this helicity amplitude squared to return all the
// same values as t_u_M_R_gqb above, TIMES a further factor tk*uk!
Energy2 VVHardGenerator::t_u_M_R_gqb_hel_amp(realVVKinematics R) const {
  using namespace ThePEG::Helicity;

//   gqb_hel_amps_.reset(ProductionMatrixElement(PDT::Spin1,PDT::Spin1Half,
// 					      PDT::Spin1,PDT::Spin1,
// 					      PDT::Spin1Half));
  
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

  VectorWaveFunction gPolarization(R.p1r(),p2data,incoming);
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
      for(unsigned int ix=0;ix<3;++ix) tc.push_back(getParticleData(1+2*ix));
    else
      for(unsigned int ix=0;ix<3;++ix) tc.push_back(getParticleData(2+2*ix));
  }
  else if(k1data->id()==23&&k2data->id()==23)      tc.push_back(p2data);
  else if(abs(k1data->id())==24&&k2data->id()==23) tc.push_back(kdata->CC());

  // Loop over helicities summing the relevant diagrams
  for(unsigned int p1hel=0;p1hel<2;++p1hel) {
    for(unsigned int p2hel=0;p2hel<2;++p2hel) {
      for(unsigned int k1hel=0;k1hel<3;++k1hel) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
	  for(unsigned int khel=0;khel<2;++khel) {
	    vector<Complex> diagrams;
	    SpinorBarWaveFunction p1_p2 = ffg->evaluate(QCDScale_,5,p2data,qbin[p2hel],g[p1hel]);
	    SpinorWaveFunction    p1_k  = ffg->evaluate(QCDScale_,5,kdata->CC(),qbout[khel],g[p1hel]);
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
// 	    gqb_hel_amps_(p1hel,p2hel,k1hel,k2hel,khel) = hel_amp;
	    sum_hel_amps_sqr += norm(hel_amp);
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
double VVHardGenerator::lo_me() const {
  using namespace ThePEG::Helicity;

  double sum_hel_amps_sqr(0.);

  tcPDPtr p1data(quark_->dataPtr());
  tcPDPtr p2data(antiquark_->dataPtr());
  tcPDPtr k1data(V1_->dataPtr());
  tcPDPtr k2data(V2_->dataPtr());
  if(k1data->id()==-24&&k2data->id()==24) swap(k1data,k2data); // Should never actually occur.

  SpinorWaveFunction qSpinor(B_.p1b(),p1data,incoming);
  SpinorBarWaveFunction qbSpinor(B_.p2b(),p2data,incoming);
  vector<SpinorWaveFunction> q;
  vector<SpinorBarWaveFunction> qb;
  for(unsigned int ix=0;ix<2;ix++) {
    qSpinor.reset(ix);
    qbSpinor.reset(ix);
    q.push_back(qSpinor);
    qb.push_back(qbSpinor);
  }

  VectorWaveFunction v1Polarization(B_.k1b(),k1data,outgoing);
  VectorWaveFunction v2Polarization(B_.k2b(),k2data,outgoing);
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
      for(unsigned int ix=0;ix<3;++ix) tc.push_back(getParticleData(1+2*ix));
    else
      for(unsigned int ix=0;ix<3;++ix) tc.push_back(getParticleData(2+2*ix));
  }
  else if(k1data->id()==23&&k2data->id()==23)      tc.push_back(p1data);
  else if(abs(k1data->id())==24&&k2data->id()==23) tc.push_back(p2data);

  // Loop over helicities summing the relevant diagrams
  for(unsigned int p1hel=0;p1hel<2;++p1hel) {
    for(unsigned int p2hel=0;p2hel<2;++p2hel) {
      for(unsigned int k1hel=0;k1hel<3;++k1hel) {
	for(unsigned int k2hel=0;k2hel<3;++k2hel) {
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
	  sum_hel_amps_sqr += norm(hel_amp);
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

// // Pseudo-code for including spin correlations:
// void VVHardGenerator::constructVertex(tSubProPtr sub) {
//   //  if(!spinCorrelations()) return;
//   SpinfoPtr spin[5];
//   // Extract the particles in the hard process
//   ParticleVector hard;
//   hard.push_back(sub->incoming().first);  // p1 
//   hard.push_back(sub->incoming().second); // p2
//   hard.push_back(sub->outgoing()[0]);     // k1
//   hard.push_back(sub->outgoing()[1]);     // k2
//   hard.push_back(sub->outgoing()[2]);     // k
//   // Get the spin info objects
//   for(unsigned int ix=0;ix<5;++ix)
//     spin[ix]=dynamic_ptr_cast<SpinfoPtr>(hard[]->spinInfo());
//   // construct the vertex
//   HardVertexPtr hardvertex=new_ptr(HardVertex());
//   // set the matrix element for the vertex
//   hardvertex->ME(qqb_hel_amps_);
//   // set the pointers and to and from the vertex
//   for(unsigned int ix=0;ix<5;++ix)
//     spin[ix]->setProductionVertex(hardvertex);
// }

