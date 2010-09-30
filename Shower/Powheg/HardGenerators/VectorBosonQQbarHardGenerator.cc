// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorBosonQQbarHardGenerator class.
//

#include "VectorBosonQQbarHardGenerator.h"
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
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/Switch.h"

using namespace Herwig;

void VectorBosonQQbarHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS << _alphaEM << _gamma << _gluon << ounit( _Ptmin, GeV );
}

void VectorBosonQQbarHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> _alphaEM >> _gamma >> _gluon >> iunit( _Ptmin, GeV );
}

ClassDescription<VectorBosonQQbarHardGenerator> 
VectorBosonQQbarHardGenerator::initVectorBosonQQbarHardGenerator;
// Definition of the static class description member.

void VectorBosonQQbarHardGenerator::Init() {

  static ClassDocumentation<VectorBosonQQbarHardGenerator> documentation
    ("The VectorBosonQQbarHardGenerator class generates the hardest emission for"
     "vector boson decays to fermion-antifermion events in the POWHEG approach");

  static Reference<VectorBosonQQbarHardGenerator,ShowerAlpha> interfaceShowerAlphaQCD
    ("ShowerAlphaQCD",
     "The object calculating the strong coupling constant",
     &VectorBosonQQbarHardGenerator::_alphaS, false, false, true, false, false);
  

  static Reference<VectorBosonQQbarHardGenerator,ShowerAlpha> interfaceShowerAlphaQED
    ("ShowerAlphaQED",
     "The electromagnetic coupling for the QED corrections",
     &VectorBosonQQbarHardGenerator::_alphaEM, false, false, true, false, false);

  static Parameter<VectorBosonQQbarHardGenerator, Energy> interfacePtMin
    ("minPt",
     "The pt cut on hardest emision generation",
     &VectorBosonQQbarHardGenerator::_Ptmin, GeV, 1.*GeV, 0*GeV, 100000.0*GeV,
     false, false, Interface::limited);
}

void VectorBosonQQbarHardGenerator::doinit() {
  HardestEmissionGenerator::doinit();
  _gluon = getParticleData(ParticleID::g);
  _gamma = getParticleData(ParticleID::gamma);
}


bool VectorBosonQQbarHardGenerator::canHandle(ShowerTreePtr tree) {
  // one incoming gauge boson
  if(tree->incomingLines().size()!=1) return false;
  // what on earth does this do??
  if((tree->incomingLines().begin()->first->id()==22)&&
     (tree->incomingLines().begin()->first->progenitor()->id()==23)) {
    generator()->log() << "testing in dodgy if\n";
    return false;
  }
  // check outgoing quarks or leptons
  map<ShowerProgenitorPtr,tShowerParticlePtr> outgoing=tree->outgoingLines();
  // must be two
  if(outgoing.size()!=2) return false;
  // particle and antiparticle
  if(outgoing.begin()->first->progenitor()->id()!=
     -1*outgoing.rbegin()->first->progenitor()->id())     return false;
  int id = abs(outgoing.begin()->first->progenitor()->id());
  // quark or charged lepton 
  if( id <= 6 || (id >= 11 && id <= 15 && id%2 == 1))  return true;
  return false;
}

HardTreePtr VectorBosonQQbarHardGenerator::generateHardest(ShowerTreePtr tree) {
  // Get the progenitors: Q and Qbar.
  ShowerProgenitorPtr 
    QProgenitor    = tree->outgoingLines().begin()->first,
    QbarProgenitor = tree->outgoingLines().rbegin()->first;
  if(QProgenitor->id()<0) swap( QProgenitor, QbarProgenitor );
  _partons.resize(2);
  _partons[0] = QProgenitor->progenitor()->dataPtr();
  _partons[1] = QbarProgenitor->progenitor()->dataPtr();
  // momentum of the partons
  _quark.resize(2);
  _quark[0]=QProgenitor->copy()->momentum();
  _quark[1]=QbarProgenitor->copy()->momentum();
  // Set the existing mass entries in partons 5 vectors with the
  // once and for all.
  _quark[0].setMass(_partons[0]->mass());
  _quark[1].setMass(_partons[1]->mass());
  _gauge.setMass(0.*MeV);
  // Get the gauge boson.
  _boson = tree->incomingLines().begin()->first->copy();
  // Get the gauge boson mass.
  _s = (_quark[0]+_quark[1]).m2();
  // Generate emission and set _quark[0,1] and _gauge to be the 
  // momenta of q, qbar and g after the hardest emission:
  if(!getEvent()) {
    QProgenitor->maximumpT(_Ptmin);
    QbarProgenitor->maximumpT(_Ptmin);
    return HardTreePtr();
  }
  // Get the list of possible branchings.
  BranchingList branchings = 
    evolver()->splittingGenerator()->finalStateBranchings();
  // PDG codes of the partons
  int q_id   (abs(QProgenitor->progenitor()->id()   ));
  int qbar_id(abs(QbarProgenitor->progenitor()->id()));
  // Find the sudakovs for the q/qbar->q/qbarg branchings.
  SudakovPtr q_sudakov,qbar_sudakov;
  // quark
  long q_index(q_id),qbar_index(qbar_id); 
  for(BranchingList::const_iterator cit = branchings.lower_bound(q_index);
      cit != branchings.upper_bound(q_index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0]==q_id&&ids[1]==q_id&&
       ((_inter==ShowerInteraction::QCD&&ids[2]==ParticleID::g    )||
	(_inter==ShowerInteraction::QED&&ids[2]==ParticleID::gamma))) {
      q_sudakov = cit->second.first;
      break; 	    
    }
  }
  // antiquark
  for(BranchingList::const_iterator cit = branchings.lower_bound(qbar_index);
      cit != branchings.upper_bound(qbar_index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0]==qbar_id&&ids[1]==qbar_id&&
       ((_inter==ShowerInteraction::QCD&&ids[2]==ParticleID::g    )||
	(_inter==ShowerInteraction::QED&&ids[2]==ParticleID::gamma))) {
      qbar_sudakov = cit->second.first;
      break; 	    
    }
  }
  // Check sudakovs got created:
  if(!q_sudakov||!qbar_sudakov) throw Exception() 
    << "No Sudakov for the hard emission in "
    << "VectorBosonQQbarHardGenerator::generateHardest()" 
    << Exception::runerror;
  // Ensure the energies are greater than the constituent masses:
  for (int i=0; i<2; i++)
    if (_quark[i].e() < _partons[i]->constituentMass()) return HardTreePtr();
  if (_inter==ShowerInteraction::QCD && 
      _gauge.e() < _gluon->constituentMass()) return HardTreePtr();
  // set masses
  _quark[0].setMass( _partons[0]->mass() );
  _quark[1].setMass( _partons[1]->mass() );
  _gauge.setMass( ZERO );
  // assign the emitter based on evolution scales
  _iemitter   = _quark[0]*_gauge>_quark[1]*_gauge ? 1 : 0;
  _ispectator = _iemitter==1              ? 0 : 1; 
  // Set the sudakov for the q/qbar->q/qbarg branching.
  SudakovPtr sudakov = _iemitter == 0 ? q_sudakov : qbar_sudakov;
  // Make the particles for the HardTree:
  ShowerParticlePtr emitter(new_ptr(ShowerParticle(_partons[_iemitter],true)));
  ShowerParticlePtr spectator(new_ptr(ShowerParticle(_partons[_ispectator],true)));
  ShowerParticlePtr gauge(new_ptr(ShowerParticle(_inter == ShowerInteraction::QCD ?
						 _gluon : _gamma ,true)));
  ShowerParticlePtr vboson(new_ptr(ShowerParticle(_boson->dataPtr(),false)));
  ShowerParticlePtr parent(new_ptr(ShowerParticle(_partons[_iemitter],true)));
  emitter->set5Momentum(_quark[_iemitter]); 
  spectator->set5Momentum(_quark[_ispectator]);  
  gauge->set5Momentum(_gauge);  
  vboson->set5Momentum(_boson->momentum());  
  Lorentz5Momentum parentMomentum(_quark[_iemitter]+_gauge);
  parentMomentum.setMass(_partons[_iemitter]->mass());
  parent->set5Momentum(parentMomentum);
  // Create the vectors of HardBranchings to create the HardTree:
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  // Incoming boson:
  spaceBranchings.push_back(new_ptr(HardBranching(vboson,SudakovPtr(),
						  HardBranchingPtr(),
						  HardBranching::Incoming)));
  // Outgoing particles from hard emission:
  HardBranchingPtr spectatorBranch(new_ptr(HardBranching(spectator,
							 SudakovPtr(),HardBranchingPtr(),
							 HardBranching::Outgoing)));
  HardBranchingPtr emitterBranch(new_ptr(HardBranching(parent,
						       sudakov,HardBranchingPtr(),
						       HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(emitter, 
						SudakovPtr(),HardBranchingPtr(),
						HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(gauge,
						SudakovPtr(),HardBranchingPtr(),
						HardBranching::Outgoing)));
  if(_iemitter==0) {
    allBranchings.push_back(emitterBranch);
    allBranchings.push_back(spectatorBranch);
  } 
  else {
    allBranchings.push_back( spectatorBranch );
    allBranchings.push_back( emitterBranch );
  }
  // Add incoming boson to allBranchings
  allBranchings.push_back( spaceBranchings.back() );
  // Make the HardTree from the HardBranching vectors.
  HardTreePtr nasontree = new_ptr(HardTree(allBranchings,spaceBranchings,
					   _inter));
  // Set the maximum pt for all other emissions
  Energy ptveto(_pt);
  QProgenitor   ->maximumpT(ptveto);
  QbarProgenitor->maximumpT(ptveto);
  // Connect the particles with the branchings in the HardTree
  nasontree->connect( QProgenitor->progenitor(), allBranchings[0] );
  nasontree->connect( QbarProgenitor->progenitor(), allBranchings[1] );
  // Create the two colour lines and connect the particles:
  ColinePtr blueLine  = new_ptr(ColourLine());
  ColinePtr greenLine = new_ptr(ColourLine());
  blueLine->addColoured(emitter,_iemitter);
  blueLine->addColoured(gauge,_ispectator);
  greenLine->addColoured(gauge,_iemitter);
  greenLine->addColoured(parent,_iemitter);
  greenLine->addColoured(spectator,_ispectator);
  ShowerParticleVector particles;
  for(set<HardBranchingPtr>::iterator cit=nasontree->branchings().begin();
      cit!=nasontree->branchings().end();++cit) {
    particles.push_back((*cit)->branchingParticle());
  }
  // set the evolution partners and scales
  evolver()->showerModel()->partnerFinder()->
    setInitialEvolutionScales(particles,true,_inter,true);
  // Calculate the shower variables:
  evolver()->showerModel()->kinematicsReconstructor()->
    deconstructDecayJets(nasontree,evolver(),ShowerInteraction::QCD);
  // Return the HardTree
  return nasontree;
}

bool VectorBosonQQbarHardGenerator::getEvent() {
  // max pT
  Energy ptmax = 0.5*sqrt(_s);
  // Define over valued y_max & y_min according to the associated pt_min cut.
  double ymax  =  acosh(ptmax/_Ptmin);
  double ymin  = -acosh(ptmax/_Ptmin);
  // loop over possible QED and QCD interactions
  Energy pT[2]={-GeV,-GeV};
  double  y[2];
  for(unsigned int interaction=0;interaction<2;++interaction) { 
    // pt of the emmission
    Energy lastpt(ptmax);
    // rapidity of the emission
    double lasty(0.);
    // pefactor
    double prefactor(0.);
    if(interaction==0&&_partons[0]->coloured()) {
      prefactor = 2.*_alphaS->overestimateValue()*(4./3.)/Constants::pi;
    }
    else if(interaction==1&&_partons[0]->charged()) {
      prefactor = 2.*_alphaEM->overestimateValue()*
	sqr(double(_partons[0]->iCharge())/3.)/Constants::pi;
    }
    // if not allowed continue
    if(prefactor==0.) continue;
    // loop to generate the pt and rapidity
    bool reject;
    do {
      // don't reject the emission
      reject = true;
      // generator pt
      lastpt *= pow(UseRandom::rnd(),1./(prefactor*(ymax-ymin)));
      if(lastpt<_Ptmin) {
	lastpt=-GeV;
	break;
      }
      lasty   = ymin + UseRandom::rnd()*(ymax-ymin);
      // calculate x1 and x2
      double x1 = 1.-lastpt/sqrt(_s)*exp( lasty);
      double x2 = 1.-lastpt/sqrt(_s)*exp(-lasty);
      if(x1<0.||x1>1.||x2<0.||x2>1.) {
	reject=true;
	continue;
      }
      // matrix element weight
      double wgt = (sqr(x1)+sqr(x2))/(1.-x1)/(1.-x2)*0.5*sqr(lastpt)/_s;
      if(interaction==0)      wgt *= _alphaS ->ratio(sqr(lastpt));
      else if(interaction==1) wgt *= _alphaEM->ratio(sqr(lastpt));
      if(wgt>1.) { 
	generator()->log() << "VectorBosonQQbarHardGenerator::getEvent() " 
			   << "excess weight " << wgt << "\n";
      }
      reject = UseRandom::rnd()>wgt;
    }
    while(reject);
    pT[interaction] = lastpt;
    y [interaction] = lasty;
  }
  if(pT[0]>pT[1]) {
    _inter = ShowerInteraction::QCD;
    _pt = pT[0];
    _y  =  y[0];
  }
  else {
    _inter = ShowerInteraction::QED;
    _pt = pT[1];
    _y  =  y[1];
  }
  // no emission
  if(_pt<ZERO) return false;
  // x values
  _xq  = 1.-_pt/sqrt(_s)*exp( _y);
  _xqb = 1.-_pt/sqrt(_s)*exp(-_y);
  _xg  = 2.-_xq-_xqb;
  // select emitter and spectator
  UseRandom::rnd()<(sqr(_xq)/(sqr(_xq)+sqr(_xqb))) ? 
    _iemitter = 1: _iemitter = 0 ; 
  _ispectator = !_iemitter;
  //construct vectors in com z frame
  return constructVectors();
}

bool VectorBosonQQbarHardGenerator::constructVectors(){
  using Constants::twopi;
  using Constants::pi;
  // Find the boost from the lab to the c.o.m with the spectator 
  // along the -z axis, and then invert it.
  LorentzRotation eventFrame( ( _quark[0] + _quark[1] ).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*_quark[_ispectator];
  eventFrame.rotateZ( -spectator.phi() );
  eventFrame.rotateY( -spectator.theta() - pi );
  eventFrame.invert();

  Energy rts = sqrt(_s);

  // Extract the reduced (constituent) masses:
  double mu_e(_quark[_iemitter].mass()/rts);
  double mu_s(_quark[_ispectator].mass()/rts);
  double mu_g(_gauge.mass()/rts);

  // Get masses to avoid floating point errors later.
  Energy init_g_mass(_gauge.mass());
  Energy init_e_mass(_quark[_iemitter].mass());
  Energy init_s_mass(_quark[_ispectator].mass());
  // Get the Dalitz variables:
  double xe   = _iemitter==0   ? _xq : _xqb;
  double xs   = _ispectator==0 ? _xq : _xqb;
  double xg   = _xg;
  double b_xe2  = (xe + 2.*mu_e )*(xe - 2.*mu_e );
  double b_xs2 =  (xs + 2.*mu_s )*(xs - 2.*mu_s );
  double b_xg2  = (xg + 2.*mu_g )*(xg - 2.*mu_g );
  if(b_xe2 <0. || b_xs2<0.|| b_xg2 <0.) return false;
  double b_xe = sqrt(b_xe2 );
  double b_xs = sqrt(b_xs2 );
  double b_xg = sqrt(b_xg2 );
  if(xe > (1.+(mu_e - mu_s - mu_g )*(mu_e + mu_s + mu_g ))) return false;
  if(xs > (1.+(mu_s - mu_e - mu_g )*(mu_s + mu_e + mu_g ))) return false;
  // Get the cosines and sines of emitter w.r.t spectator:
  double c_se,s_se;
  c_se  =  0.5*(b_xe2+b_xs2-b_xg2)/(b_xs*b_xe);
  s_se  = sqrt(1.-sqr(c_se));

  if(abs(c_se)>1.||abs(s_se)>1.) return false;

  // momentum of emitter
  _quark[_iemitter].setT( 0.5*rts * xe);
  _quark[_iemitter].setX( 0.5*rts * b_xe*s_se*cos(_phi));
  _quark[_iemitter].setY( 0.5*rts * b_xe*s_se*sin(_phi));
  _quark[_iemitter].setZ( 0.5*rts * b_xe*c_se);

  //spectator
  _quark[_ispectator].setT( 0.5*rts * xs);
  _quark[_ispectator].setX( 0.*MeV);
  _quark[_ispectator].setY( 0.*MeV);
  _quark[_ispectator].setZ(-0.5*rts * b_xs);

  // momentum of gluon
  _gauge=-_quark[0]-_quark[1];
  _gauge.setT(sqrt(_s)+_gauge.t());
  _gauge.setMass(init_g_mass);
  _gauge.rescaleRho();

  // boost constructed vectors into the event frame
  _quark[0] = eventFrame * _quark[0];
  _quark[1] = eventFrame * _quark[1];
  _gauge        = eventFrame * _gauge;
  // need to reset masses because for whatever reason the boost  
  // touches the mass component of the five-vector and can make  
  // zero mass objects acquire a floating point negative mass(!).
  _gauge.setMass(init_g_mass);
  _quark[_iemitter].setMass(init_e_mass);
  _quark[_ispectator].setMass(init_s_mass);
  _gauge.rescaleRho();
  _quark[_iemitter].rescaleRho();
  _quark[_ispectator].rescaleRho();
  return true;
}
