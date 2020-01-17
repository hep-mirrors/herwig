// -*- C++ -*-
//
// KinematicsReconstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KinematicsReconstructor class.
//

#include "KinematicsReconstructor.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "Herwig/Shower/QTilde/Base/PartnerFinder.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Shower/QTilde/SplittingFunctions/SplittingFunction.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/EventRecord/ColourLine.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Shower/QTilde/QTildeShowerHandler.h"
#include <cassert>
#include "KinematicsReconstructor.tcc"

using namespace Herwig;

DescribeClass<KinematicsReconstructor,Interfaced>
describeKinematicsReconstructor("Herwig::KinematicsReconstructor", "HwShower.so");

namespace {

/**
 *  Struct to order the jets in off-shellness
 */
struct JetOrdering {

  bool operator() (const JetKinStruct & j1, const JetKinStruct & j2) {
    Energy diff1 = j1.q.m()-j1.p.m();
    Energy diff2 = j2.q.m()-j2.p.m();
    if(diff1!=diff2) {
      return diff1>diff2;
    }
    else if( j1.q.e() != j2.q.e() )
      return j1.q.e()>j2.q.e();
    else
      return j1.parent->uniqueId>j2.parent->uniqueId;
  }
};

}

void KinematicsReconstructor::persistentOutput(PersistentOStream & os) const {
  os << _reconopt << _initialBoost << ounit(_minQ,GeV) << _noRescale 
     << _noRescaleVector << _initialStateReconOption << _finalFinalWeight;
}

void KinematicsReconstructor::persistentInput(PersistentIStream & is, int) {
  is >> _reconopt >> _initialBoost >> iunit(_minQ,GeV) >> _noRescale 
     >> _noRescaleVector >> _initialStateReconOption >> _finalFinalWeight;
}

void KinematicsReconstructor::Init() {

  static ClassDocumentation<KinematicsReconstructor> documentation
    ( "This class is responsible for the kinematics reconstruction of the showering,",
      " including the kinematics reshuffling necessary to compensate for the recoil"
      "of the emissions." );

  static Switch<KinematicsReconstructor,unsigned int> interfaceReconstructionOption
    ("ReconstructionOption",
     "Option for the kinematics reconstruction",
     &KinematicsReconstructor::_reconopt, 0, false, false);
  static SwitchOption interfaceReconstructionOptionGeneral
    (interfaceReconstructionOption,
     "General",
     "Use the general solution which ignores the colour structure for all processes",
     0);
  static SwitchOption interfaceReconstructionOptionColour
    (interfaceReconstructionOption,
     "Colour",
     "Use the colour structure of the process to determine the reconstruction procedure.",
     1);
  static SwitchOption interfaceReconstructionOptionColour2
    (interfaceReconstructionOption,
     "Colour2",
     "Make the most use possible of the colour structure of the process to determine the reconstruction procedure. "
     "Start with FF, then IF then II colour connections",
     2);
  static SwitchOption interfaceReconstructionOptionColour3
    (interfaceReconstructionOption,
     "Colour3",
     "Make the most use possible of the colour structure of the process to determine the reconstruction procedure. "
     "Do the colour connections in order of the pT's emitted in the shower starting with the hardest."
     " The colour partner is fully reconstructed at the same time.",
     3);
  static SwitchOption interfaceReconstructionOptionColour4
    (interfaceReconstructionOption,
     "Colour4",
     "Make the most use possible of the colour structure of the process to determine the reconstruction procedure. "
     "Do the colour connections in order of the pT's emitted in the shower starting with the hardest, while leaving"
     " the colour partner on mass-shell",
     4);

  static Parameter<KinematicsReconstructor,Energy> interfaceMinimumQ2
    ("MinimumQ2",
     "The minimum Q2 for the reconstruction of initial-final systems",
     &KinematicsReconstructor::_minQ, GeV, 0.001*GeV, 1e-6*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static RefVector<KinematicsReconstructor,ParticleData> interfaceNoRescale
    ("NoRescale",
     "Particles which shouldn't be rescaled to be on shell by the shower",
     &KinematicsReconstructor::_noRescaleVector, -1, false, false, true, false, false);

  static Switch<KinematicsReconstructor,unsigned int> interfaceInitialInitialBoostOption
    ("InitialInitialBoostOption",
     "Option for how the boost from the system before ISR to that after ISR is applied.",
     &KinematicsReconstructor::_initialBoost, 0, false, false);
  static SwitchOption interfaceInitialInitialBoostOptionOneBoost
    (interfaceInitialInitialBoostOption,
     "OneBoost",
     "Apply one boost from old CMS to new CMS",
     0);
  static SwitchOption interfaceInitialInitialBoostOptionLongTransBoost
    (interfaceInitialInitialBoostOption,
     "LongTransBoost",
     "First apply a longitudinal and then a transverse boost",
     1);

  static Deleted<KinematicsReconstructor> delFinalStateReconOption
    ("FinalStateReconOption", "The old default (0) is now the only choice");
  
  static Switch<KinematicsReconstructor,unsigned int> interfaceInitialStateReconOption
    ("InitialStateReconOption",
     "Option for the reconstruction of initial state radiation",
     &KinematicsReconstructor::_initialStateReconOption, 0, false, false);
  static SwitchOption interfaceInitialStateReconOptionRapidity
    (interfaceInitialStateReconOption,
     "Rapidity",
     "Preserve shat and rapidity",
     0);
  static SwitchOption interfaceInitialStateReconOptionLongitudinal
    (interfaceInitialStateReconOption,
     "Longitudinal",
     "Preserve longitudinal momentum",
     1);
  static SwitchOption interfaceInitialStateReconOptionSofterFraction
    (interfaceInitialStateReconOption,
     "SofterFraction",
     "Preserve the momentum fraction of the parton which has emitted softer.",
     2);

  static Switch<KinematicsReconstructor,bool> interfaceFinalFinalWeight
    ("FinalFinalWeight",
     "Apply kinematic rejection weight for final-states",
     &KinematicsReconstructor::_finalFinalWeight, false, false, false);
  static SwitchOption interfaceFinalFinalWeightNo
    (interfaceFinalFinalWeight,
     "No",
     "Don't apply the weight",
     false);
  static SwitchOption interfaceFinalFinalWeightYes
    (interfaceFinalFinalWeight,
     "Yes",
     "Apply the weight",
     true);

}

void KinematicsReconstructor::doinit() {
  Interfaced::doinit();
  _noRescale = set<cPDPtr>(_noRescaleVector.begin(),_noRescaleVector.end());
}

bool KinematicsReconstructor::
reconstructTimeLikeJet(const tShowerParticlePtr particleJetParent) const {
  assert(particleJetParent);
  bool emitted=true;
  // if this is not a fixed point in the reconstruction
  if( !particleJetParent->children().empty() ) {
    // if not a reconstruction fixpoint, dig deeper for all children:
    for ( ParticleVector::const_iterator cit = 
	    particleJetParent->children().begin();
	  cit != particleJetParent->children().end(); ++cit )
      reconstructTimeLikeJet(dynamic_ptr_cast<ShowerParticlePtr>(*cit));
  }
  // it is a reconstruction fixpoint, ie kinematical data has to be available 
  else {
    // check if the parent was part of the shower
    ShowerParticlePtr jetGrandParent;
    if(!particleJetParent->parents().empty())
      jetGrandParent= dynamic_ptr_cast<ShowerParticlePtr>
	(particleJetParent->parents()[0]);
    // update if so
    if (jetGrandParent) {
      if (jetGrandParent->showerKinematics()) {
	if(particleJetParent->id()==_progenitor->id()&&
	   !_progenitor->data().stable()&&abs(_progenitor->data().id())!=ParticleID::tauminus) {
	  jetGrandParent->showerKinematics()->reconstructLast(particleJetParent,
							      _progenitor->mass());
	}
	else {
	  jetGrandParent->showerKinematics()->reconstructLast(particleJetParent);
	}
      }
    }
    // otherwise
    else {
      Energy dm = ShowerHandler::currentHandler()->retConstituentMasses()?
	particleJetParent->data().constituentMass():
	particleJetParent->data().mass();
      if (abs(dm-particleJetParent->momentum().m())>0.001*MeV
	  &&(particleJetParent->dataPtr()->stable() || abs(particleJetParent->id())==ParticleID::tauminus)
	  &&particleJetParent->id()!=ParticleID::gamma
	  &&_noRescale.find(particleJetParent->dataPtr())==_noRescale.end()) {
	Lorentz5Momentum dum =  particleJetParent->momentum();
	dum.setMass(dm);
	dum.rescaleEnergy();
	if(abs(particleJetParent->id())==15&&particleJetParent->spinInfo()) {
	  if(particleJetParent->spinInfo()->isNear(particleJetParent->momentum())) {
	    particleJetParent->spinInfo()->SpinInfo::transform(dum,LorentzRotation());
	  }
	}
	particleJetParent->set5Momentum(dum);
      } 
      else {
	emitted=false;
      }
    }
  }
  // recursion has reached an endpoint once, ie we can reconstruct the
  // kinematics from the children.
  if( !particleJetParent->children().empty() ) 
    particleJetParent->showerKinematics()
      ->reconstructParent( particleJetParent, particleJetParent->children() );
  return emitted;
}

bool KinematicsReconstructor::
reconstructHardJets(ShowerTreePtr hard,
		    const map<tShowerProgenitorPtr,
		    pair<Energy,double> > & intrinsic,
		    ShowerInteraction type,
		    bool switchRecon) const {
  _currentTree = hard;
  _intrinsic=intrinsic;
  // extract the particles from the ShowerTree
  vector<ShowerProgenitorPtr> ShowerHardJets=hard->extractProgenitors();
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    _boosts[ShowerHardJets[ix]->progenitor()] = vector<LorentzRotation>();
  }
  for(map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator
	tit  = _currentTree->treelinks().begin();
      tit != _currentTree->treelinks().end();++tit) {
    _treeBoosts[tit->first] = vector<LorentzRotation>();
  }
  try {
    // old recon method, using new member functions
    if(_reconopt == 0 || switchRecon ) {
      reconstructGeneralSystem(ShowerHardJets);
    }
    // reconstruction based on coloured systems
    else if( _reconopt == 1) {
      reconstructColourSinglets(ShowerHardJets,type);
    }
    // reconstruction of FF, then IF, then II
    else if( _reconopt == 2) {
      reconstructFinalFirst(ShowerHardJets);
    }
    // reconstruction based on coloured systems
    else if( _reconopt == 3 || _reconopt == 4) {
      reconstructColourPartner(ShowerHardJets);
    }
    else
      assert(false);
  }
  catch(KinematicsReconstructionVeto) {
    _progenitor=tShowerParticlePtr();
    _intrinsic.clear();
    for(map<tPPtr,vector<LorentzRotation> >::const_iterator bit=_boosts.begin();bit!=_boosts.end();++bit) {
      for(vector<LorentzRotation>::const_reverse_iterator rit=bit->second.rbegin();rit!=bit->second.rend();++rit) {
	LorentzRotation rot = rit->inverse();
	bit->first->transform(rot);
      }
    }
    _boosts.clear();
    for(map<tShowerTreePtr,vector<LorentzRotation> >::const_iterator bit=_treeBoosts.begin();bit!=_treeBoosts.end();++bit) {
      for(vector<LorentzRotation>::const_reverse_iterator rit=bit->second.rbegin();rit!=bit->second.rend();++rit) {
	LorentzRotation rot = rit->inverse();
	bit->first->transform(rot,false);
      }
    }
    _currentTree = tShowerTreePtr();
    _treeBoosts.clear();
    return false;
  }
  catch (Exception & ex) {
    _progenitor=tShowerParticlePtr();
    _intrinsic.clear();
    _currentTree = tShowerTreePtr();
    _boosts.clear();
    _treeBoosts.clear();
    throw ex;
  }
  _progenitor=tShowerParticlePtr();
  _intrinsic.clear();
  // ensure x<1
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=hard->incomingLines().begin();cit!=hard->incomingLines().end();++cit) {
    tPPtr parent = cit->first->progenitor();
    while (!parent->parents().empty()) {
      parent = parent->parents()[0];
    }
    tPPtr hadron;
    if ( cit->first->original()->parents().empty() ) {
      hadron = cit->first->original();
    } 
    else {
      hadron = cit->first->original()->parents()[0];
    }
    if( ! (hadron->id() == parent->id() && hadron->children().size() <= 1)
       && parent->momentum().rho() > hadron->momentum().rho()) {
      _progenitor=tShowerParticlePtr();
      _intrinsic.clear();
      for(map<tPPtr,vector<LorentzRotation> >::const_iterator bit=_boosts.begin();bit!=_boosts.end();++bit) {
	for(vector<LorentzRotation>::const_reverse_iterator rit=bit->second.rbegin();rit!=bit->second.rend();++rit) {
	  LorentzRotation rot = rit->inverse();
	  bit->first->transform(rot);
	}
      }
      _boosts.clear();
      for(map<tShowerTreePtr,vector<LorentzRotation> >::const_iterator bit=_treeBoosts.begin();bit!=_treeBoosts.end();++bit) {
	for(vector<LorentzRotation>::const_reverse_iterator rit=bit->second.rbegin();rit!=bit->second.rend();++rit) {
	  LorentzRotation rot = rit->inverse();
	  bit->first->transform(rot,false);
	}
      }
      _currentTree = tShowerTreePtr();
      _treeBoosts.clear();
      return false;
    }
  }
  _boosts.clear();
  _treeBoosts.clear();
  _currentTree = tShowerTreePtr();
  return true;
}

double 
KinematicsReconstructor::solveKfactor(const Energy & root_s, 
				      const JetKinVect & jets) const {
  Energy2 s = sqr(root_s);
  // must be at least two jets
  if ( jets.size() < 2) throw KinematicsReconstructionVeto();
  // sum of jet masses must be less than roots
  if(momConsEq( 0.0, root_s, jets )>ZERO) throw KinematicsReconstructionVeto();
  // if two jets simple solution
  if ( jets.size() == 2 ) {
    static const Energy2 eps = 1.0e-4 * MeV2;
    if ( sqr(jets[0].p.x()+jets[1].p.x()) < eps &&
	 sqr(jets[0].p.y()+jets[1].p.y()) < eps &&
	 sqr(jets[0].p.z()+jets[1].p.z()) < eps ) {
      Energy test = (jets[0].p+jets[1].p).vect().mag();
      if(test > 1.0e-4 * MeV) throw KinematicsReconstructionVeto();
      if ( jets[0].p.vect().mag2() < eps ) throw KinematicsReconstructionVeto();
      Energy2 m1sq(jets[0].q.m2()),m2sq(jets[1].q.m2());
      return sqrt( ( sqr(s - m1sq - m2sq) - 4.*m1sq*m2sq )
		   /(4.*s*jets[0].p.vect().mag2()) );
    } 
    else throw KinematicsReconstructionVeto();
  }
  // i.e. jets.size() > 2, numerically
  // check convergence, if it's a problem maybe use Newton iteration?
  else {
    double k1 = 0.,k2 = 1.,k = 0.; 
    if ( momConsEq( k1, root_s, jets ) < ZERO ) {
      while ( momConsEq( k2, root_s, jets ) < ZERO ) {
	k1 = k2; 
	k2 *= 2;       
      }
      while ( fabs( (k1 - k2)/(k1 + k2) ) > 1.e-10 ) { 
	if( momConsEq( k2, root_s, jets ) == ZERO ) {
	  return k2; 
	} else {
	  k = (k1+k2)/2.;
	  if ( momConsEq( k, root_s, jets ) > ZERO ) {
	    k2 = k;
	  } else {
	    k1 = k; 
	  } 
	}
      }
      return k1; 	  
    } else throw KinematicsReconstructionVeto();
  }
  throw KinematicsReconstructionVeto(); 
}

bool KinematicsReconstructor::
reconstructSpaceLikeJet( const tShowerParticlePtr p) const {
  bool emitted = true;
  tShowerParticlePtr child;
  tShowerParticlePtr parent;
  if(!p->parents().empty())
    parent = dynamic_ptr_cast<ShowerParticlePtr>(p->parents()[0]);
  if(parent) {
    emitted=true;
    reconstructSpaceLikeJet(parent);
  }
  // if branching reconstruct time-like child
  if(p->children().size()==2)
    child = dynamic_ptr_cast<ShowerParticlePtr>(p->children()[1]);
  if(p->perturbative()==0 && child) {
    dynamic_ptr_cast<ShowerParticlePtr>(p->children()[0])->
      showerKinematics()->reconstructParent(p,p->children());
    if(!child->children().empty()) {
      _progenitor=child;
      reconstructTimeLikeJet(child);
      // calculate the momentum of the particle
      Lorentz5Momentum pnew=p->momentum()-child->momentum();
      pnew.rescaleMass();
      p->children()[0]->set5Momentum(pnew);
    }
  }
  return emitted;
}

Boost KinematicsReconstructor::
solveBoostBeta( const double k, const Lorentz5Momentum & newq,
		const Lorentz5Momentum & oldp ) {
  // try something different, purely numerical first: 
  // a) boost to rest frame of newq, b) boost with kp/E
  Energy q = newq.vect().mag(); 
  Energy2 qs = sqr(q); 
  Energy2 Q2 = newq.m2(); 
  Energy kp = k*(oldp.vect().mag()); 
  Energy2 kps = sqr(kp); 

  // usually we take the minus sign, since this boost will be smaller.
  // we only require |k \vec p| = |\vec q'| which leaves the sign of
  // the boost open but the 'minus' solution gives a smaller boost
  // parameter, i.e. the result should be closest to the previous
  // result. this is to be changed if we would get many momentum
  // conservation violations at the end of the shower from a hard
  // process.
  double betam = (q*sqrt(qs + Q2) - kp*sqrt(kps + Q2))/(kps + qs + Q2);
  // move directly to 'return' 
  Boost beta = -betam*(k/kp)*oldp.vect();
  // note that (k/kp)*oldp.vect() = oldp.vect()/oldp.vect().mag() but cheaper. 
  // leave this out if it's running properly! 
  if ( betam >= 0 ) return beta;
  else              return Boost(0., 0., 0.); 
}

bool KinematicsReconstructor::
reconstructDecayJets(ShowerTreePtr decay,
		     ShowerInteraction) const {
  _currentTree = decay;
  // extract the particles from the ShowerTree
  vector<ShowerProgenitorPtr> ShowerHardJets=decay->extractProgenitors();
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    _boosts[ShowerHardJets[ix]->progenitor()] = vector<LorentzRotation>();
  }
  for(map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator
	tit  = _currentTree->treelinks().begin();
      tit != _currentTree->treelinks().end();++tit) {
    _treeBoosts[tit->first] = vector<LorentzRotation>();
  }
  try {
    bool radiated[2]={false,false};
    // find the decaying particle and check if particles radiated
    ShowerProgenitorPtr initial;
    for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
      // only consider initial-state jets
      if(ShowerHardJets[ix]->progenitor()->isFinalState()) {
	radiated[1] |=ShowerHardJets[ix]->hasEmitted();
      }
      else {
	initial=ShowerHardJets[ix];
	radiated[0]|=ShowerHardJets[ix]->hasEmitted();
      }
    }
    // find boost to the rest frame if needed
    Boost boosttorest=-initial->progenitor()->momentum().boostVector();
    double gammarest =
      initial->progenitor()->momentum().e()/
      initial->progenitor()->momentum().mass();
    // check if need to boost to rest frame
    bool gottaBoost = (boosttorest.mag() > 1e-12);
    // if initial state radiation reconstruct the jet and set up the basis vectors
    Lorentz5Momentum pjet;
    Lorentz5Momentum nvect;
    // find the partner
    ShowerParticlePtr partner = initial->progenitor()->partner();
    Lorentz5Momentum ppartner[2];
    if(partner) ppartner[0]=partner->momentum();
    // get the n reference vector
    if(partner) {
      if(initial->progenitor()->showerKinematics()) {
	nvect = initial->progenitor()->showerBasis()->getBasis()[1];
      }
      else {
	Lorentz5Momentum ppartner=initial->progenitor()->partner()->momentum();
	if(gottaBoost) ppartner.boost(boosttorest,gammarest);
	nvect = Lorentz5Momentum( ZERO,0.5*initial->progenitor()->mass()*
				  ppartner.vect().unit()); 
	nvect.boost(-boosttorest,gammarest);
      }
    }
    // if ISR
    if(radiated[0]) {
      // reconstruct the decay jet
      reconstructDecayJet(initial->progenitor());
      // momentum of decaying particle after ISR
      pjet=initial->progenitor()->momentum()
	-decay->incomingLines().begin()->second->momentum();
      pjet.rescaleMass();
    }
    // boost initial state jet and basis vector if needed
    if(gottaBoost) {
      pjet.boost(boosttorest,gammarest);
      nvect.boost(boosttorest,gammarest);
      ppartner[0].boost(boosttorest,gammarest);
    }
    // loop over the final-state particles and do the reconstruction
    JetKinVect possiblepartners;
    JetKinVect jetKinematics;
    bool atLeastOnce = radiated[0];
    LorentzRotation restboost(boosttorest,gammarest);
    Energy inmass(ZERO);
    for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
      // only consider final-state jets
      if(!ShowerHardJets[ix]->progenitor()->isFinalState()) {
	inmass=ShowerHardJets[ix]->progenitor()->mass();
	continue;
      }
      // do the reconstruction
      JetKinStruct tempJetKin;      
      tempJetKin.parent = ShowerHardJets[ix]->progenitor();
      if(ShowerHardJets.size()==2) {
	Lorentz5Momentum dum=ShowerHardJets[ix]->progenitor()->momentum();
	dum.setMass(inmass);
	dum.rescaleRho();
	tempJetKin.parent->set5Momentum(dum);
      }
      tempJetKin.p = ShowerHardJets[ix]->progenitor()->momentum();
      if(gottaBoost) tempJetKin.p.boost(boosttorest,gammarest);
      _progenitor=tempJetKin.parent;
      if(ShowerHardJets[ix]->reconstructed()==ShowerProgenitor::notReconstructed) {
	atLeastOnce |= reconstructTimeLikeJet(tempJetKin.parent);
	ShowerHardJets[ix]->reconstructed(ShowerProgenitor::done);
      }
      if(gottaBoost) deepTransform(tempJetKin.parent,restboost);
      tempJetKin.q = ShowerHardJets[ix]->progenitor()->momentum();
      jetKinematics.push_back(tempJetKin);
    }
    if(partner) ppartner[1]=partner->momentum();
    // calculate the rescaling parameters
    double k1,k2;
    Lorentz5Momentum qt;
    if(!solveDecayKFactor(initial->progenitor()->mass(),nvect,pjet,
			  jetKinematics,partner,ppartner,k1,k2,qt)) {
      for(map<tPPtr,vector<LorentzRotation> >::const_iterator bit=_boosts.begin();bit!=_boosts.end();++bit) {
	for(vector<LorentzRotation>::const_reverse_iterator rit=bit->second.rbegin();rit!=bit->second.rend();++rit) {
	  LorentzRotation rot = rit->inverse();
	  bit->first->transform(rot);
	}
      }
      _boosts.clear();
      for(map<tShowerTreePtr,vector<LorentzRotation> >::const_iterator bit=_treeBoosts.begin();bit!=_treeBoosts.end();++bit) {
	for(vector<LorentzRotation>::const_reverse_iterator rit=bit->second.rbegin();rit!=bit->second.rend();++rit) {
	  LorentzRotation rot = rit->inverse();
	  bit->first->transform(rot,false);
	}
      }
      _treeBoosts.clear();
      _currentTree = tShowerTreePtr();
      return false;
    }
    // apply boosts and rescalings to final-state jets
    for(JetKinVect::iterator it = jetKinematics.begin(); 
	it != jetKinematics.end(); ++it) {
      LorentzRotation Trafo = LorentzRotation(); 
      if(it->parent!=partner) {
	// boost for rescaling
	if(atLeastOnce) {
	  map<tShowerTreePtr,pair<tShowerProgenitorPtr,
	    tShowerParticlePtr> >::const_iterator tit;
	  for(tit  = _currentTree->treelinks().begin();
	      tit != _currentTree->treelinks().end();++tit) {
	    if(tit->second.first && tit->second.second==it->parent)
	      break;
	  }
	  if(it->parent->children().empty()&&!it->parent->spinInfo() &&
	     tit==_currentTree->treelinks().end()) {
	    Lorentz5Momentum pnew(k2*it->p.vect(),
				  sqrt(sqr(k2*it->p.vect().mag())+it->q.mass2()),
				  it->q.mass());
	    it->parent->set5Momentum(pnew);
	  }
	  else {
	    // rescaling boost can't ever work in this case
	    if(k2<0. && it->q.mass()==ZERO)
	      throw KinematicsReconstructionVeto();
	    Trafo = solveBoost(k2, it->q, it->p);
	  }
	}
	if(gottaBoost)  Trafo.boost(-boosttorest,gammarest);
	if(atLeastOnce || gottaBoost) deepTransform(it->parent,Trafo);
      }
      else {
	Lorentz5Momentum pnew=ppartner[0];
	pnew *=k1;
	pnew-=qt;
	pnew.setMass(ppartner[1].mass());
	pnew.rescaleEnergy();
	LorentzRotation Trafo=solveBoost(1.,ppartner[1],pnew);
	if(gottaBoost) Trafo.boost(-boosttorest,gammarest);
	deepTransform(partner,Trafo);
      }
    }
  }
  catch(KinematicsReconstructionVeto) {
    for(map<tPPtr,vector<LorentzRotation> >::const_iterator bit=_boosts.begin();bit!=_boosts.end();++bit) {
      for(vector<LorentzRotation>::const_reverse_iterator rit=bit->second.rbegin();rit!=bit->second.rend();++rit) {
	LorentzRotation rot = rit->inverse();
	bit->first->transform(rot);
      }
    }
    _boosts.clear();
    for(map<tShowerTreePtr,vector<LorentzRotation> >::const_iterator bit=_treeBoosts.begin();bit!=_treeBoosts.end();++bit) {
      for(vector<LorentzRotation>::const_reverse_iterator rit=bit->second.rbegin();rit!=bit->second.rend();++rit) {
	LorentzRotation rot = rit->inverse();
	bit->first->transform(rot,false);
      }
    }
    _treeBoosts.clear();
    _currentTree = tShowerTreePtr();
    return false;
  }
  catch (Exception & ex) {
    _currentTree = tShowerTreePtr();
    _boosts.clear();
    _treeBoosts.clear();
    throw ex;
  }
  _boosts.clear();
  _treeBoosts.clear();
  _currentTree = tShowerTreePtr();
  return true;
}

bool KinematicsReconstructor::
reconstructDecayJet( const tShowerParticlePtr p) const {
  if(p->children().empty()) return false;
  tShowerParticlePtr child;
  // if branching reconstruct time-like child
  child = dynamic_ptr_cast<ShowerParticlePtr>(p->children()[1]);
  if(child) {
    _progenitor=child;
    reconstructTimeLikeJet(child);
    // calculate the momentum of the particle
    Lorentz5Momentum pnew=p->momentum()-child->momentum();
    pnew.rescaleMass();
    p->children()[0]->set5Momentum(pnew);
    child=dynamic_ptr_cast<ShowerParticlePtr>(p->children()[0]);
    reconstructDecayJet(child);
    return true;
  }
  return false;
}

bool KinematicsReconstructor::
solveDecayKFactor(Energy mb, 
		  const Lorentz5Momentum & n, 
		  const Lorentz5Momentum & pjet, 
		  const JetKinVect & jetKinematics, 
		  ShowerParticlePtr partner, 
		  Lorentz5Momentum ppartner[2],
		  double & k1, double & k2,
		  Lorentz5Momentum & qt) const {
  Energy2 pjn  = partner ? pjet.vect()*n.vect()        : ZERO;
  Energy2 pcn  = partner ? ppartner[0].vect()*n.vect() : 1.*MeV2;
  Energy2 nmag = n.vect().mag2();
  Lorentz5Momentum pn = partner ? (pjn/nmag)*n : Lorentz5Momentum();
  qt=pjet-pn; qt.setE(ZERO);
  Energy2 pt2=qt.vect().mag2();
  Energy  Ejet = pjet.e();
  // magnitudes of the momenta for fast access
  vector<Energy2> pmag;
  Energy total(Ejet);
  for(unsigned int ix=0;ix<jetKinematics.size();++ix) {
    pmag.push_back(jetKinematics[ix].p.vect().mag2());
    total+=jetKinematics[ix].q.mass();
  }
  // return if no possible solution
  if(total>mb) return false;
  Energy2 pcmag=ppartner[0].vect().mag2();
  // used newton-raphson to get the rescaling
  static const Energy eps=1e-8*GeV;
  long double d1(1.),d2(1.);
  Energy roots, ea, ec, ds;
  unsigned int ix=0;
  do {
    ++ix;
    d2    = d1 + pjn/pcn;
    roots = Ejet;
    ds    = ZERO;
    for(unsigned int iy=0;iy<jetKinematics.size();++iy) {
      if(jetKinematics[iy].parent==partner) continue;
      ea = sqrt(sqr(d2)*pmag[iy]+jetKinematics[iy].q.mass2());
      roots += ea;
      ds    += d2/ea*pmag[iy];
    }
    if(partner) {
      ec = sqrt(sqr(d1)*pcmag + pt2 + ppartner[1].mass2());
      roots += ec;
      ds    += d1/ec*pcmag;
    }
    d1    += (mb-roots)/ds;
    d2     = d1 + pjn/pcn;
  }
  while(abs(mb-roots)>eps && ix<100);
  k1=d1;
  k2=d2;
  // return true if N-R succeed, otherwise false
  return ix<100;
}

bool KinematicsReconstructor::
deconstructDecayJets(HardTreePtr decay,ShowerInteraction) const {
  // extract the momenta of the particles
  vector<Lorentz5Momentum> pin;
  vector<Lorentz5Momentum> pout;
  // on-shell masses of the decay products
  vector<Energy> mon;
  Energy mbar(-GeV);
  // the hard branchings of the particles
  set<HardBranchingPtr>::iterator cit;
  set<HardBranchingPtr> branchings=decay->branchings();
  // properties of the incoming particle
  bool ISR = false;
  HardBranchingPtr initial;
  Lorentz5Momentum qisr;
  // find the incoming particle, both before and after
  // any ISR
  for(cit=branchings.begin();cit!=branchings.end();++cit){
    if((*cit)->status()==HardBranching::Incoming||
       (*cit)->status()==HardBranching::Decay) {
      // search back up isr if needed
      HardBranchingPtr branch = *cit;
      while(branch->parent()) branch=branch->parent();
      initial=branch;
      // momentum or original parent
      pin.push_back(branch->branchingParticle()->momentum());
      // ISR?
      ISR = !branch->branchingParticle()->children().empty();
      // ISR momentum
      qisr = pin.back()-(**cit).branchingParticle()->momentum();
      qisr.rescaleMass();
    }
  }
  assert(pin.size()==1);
  // compute boost to rest frame
  Boost boostv=-pin[0].boostVector();
  // partner for ISR
  ShowerParticlePtr partner;
  Lorentz5Momentum  ppartner;
  if(initial->branchingParticle()->partner()) {
    partner=initial->branchingParticle()->partner();
    ppartner=partner->momentum();
  }
  // momentum of the decay products
  for(cit=branchings.begin();cit!=branchings.end();++cit) {
    if((*cit)->status()!=HardBranching::Outgoing) continue;
    // find the mass of the particle
    // including special treatment for off-shell resonances
    // to preserve off-shell mass
    Energy mass;
    if(!(**cit).branchingParticle()->dataPtr()->stable()) {
      HardBranchingPtr branch=*cit;
      while(!branch->children().empty()) {
	for(unsigned int ix=0;ix<branch->children().size();++ix) {
	  if(branch->children()[ix]->branchingParticle()->id()==
	     (**cit).branchingParticle()->id()) {
	    branch = branch->children()[ix];
	    continue;
	  }
	}
      };
      mass = branch->branchingParticle()->mass();
    }
    else {
      mass = (**cit).branchingParticle()->dataPtr()->mass();
    }
    // if not evolution partner of decaying particle
    if((*cit)->branchingParticle()!=partner) {
      pout.push_back((*cit)->branchingParticle()->momentum());
      mon.push_back(mass);
    }
    // evolution partner of decaying particle
    else {
      mbar = mass;
    }
  }
  // boost all the momenta to the rest frame of the decaying particle
  for(unsigned int ix=0;ix<pout.size();++ix) pout[ix].boost(boostv);
  if(initial->branchingParticle()->partner()) {
    ppartner.boost(boostv);
    qisr.boost(boostv);
  }
  // compute the rescaling factors
  double k1,k2;
  if(!ISR) {
    if(partner) {
      pout.push_back(ppartner);
      mon.push_back(mbar);
    }
    k1=k2=inverseRescalingFactor(pout,mon,pin[0].mass());
    if(partner) {
      pout.pop_back();
      mon.pop_back();
    }
  }
  else {
    if(!inverseDecayRescalingFactor(pout,mon,pin[0].mass(),
				    ppartner,mbar,k1,k2)) return false;
  }
  // now calculate the p reference vectors 
  unsigned int ifinal=0;
  for(cit=branchings.begin();cit!=branchings.end();++cit) {
    if((**cit).status()!=HardBranching::Outgoing) continue;
    // for partners other than colour partner of decaying particle
    if((*cit)->branchingParticle()!=partner) {
      Lorentz5Momentum pvect = (*cit)->branchingParticle()->momentum();
      pvect.boost(boostv);
      pvect /= k1;
      pvect.setMass(mon[ifinal]);
      ++ifinal;
      pvect.rescaleEnergy();
      pvect.boost(-boostv);
      (*cit)->pVector(pvect);
      (*cit)->showerMomentum(pvect);
    }
    // for colour partner of decaying particle
    else {
      Lorentz5Momentum pvect = (*cit)->branchingParticle()->momentum();
      pvect.boost(boostv);
      Lorentz5Momentum qtotal;
      for(unsigned int ix=0;ix<pout.size();++ix) qtotal+=pout[ix];
      Lorentz5Momentum qperp = 
	qisr-(qisr.vect()*qtotal.vect())/(qtotal.vect().mag2())*qtotal;
      pvect +=qperp;
      pvect /=k2;
      pvect.setMass(mbar);
      pvect.rescaleEnergy();
      pvect.boost(-boostv);
      (*cit)->pVector(pvect);
      (*cit)->showerMomentum(pvect);
    }
  }
  // For initial-state if needed
  if(initial) {
    tShowerParticlePtr newPartner=initial->branchingParticle()->partner();
    if(newPartner) {
      tHardBranchingPtr branch;
      for( set<HardBranchingPtr>::iterator clt = branchings.begin();
	   clt != branchings.end(); ++clt ) {
	if((**clt).branchingParticle()==newPartner) {
	  initial->colourPartner(*clt);
	  branch=*clt;
	  break;
	}
      }
      Lorentz5Momentum pvect = initial->branchingParticle()->momentum();
      initial->pVector(pvect);
      Lorentz5Momentum ptemp = branch->pVector();
      ptemp.boost(boostv);
      Lorentz5Momentum nvect = Lorentz5Momentum( ZERO,
						 0.5*initial->branchingParticle()->mass()*
						 ptemp.vect().unit());
      nvect.boost(-boostv);
      initial->nVector(nvect);
    }
  }
  // calculate the reference vectors, then for outgoing particles
  for(cit=branchings.begin();cit!=branchings.end();++cit){
    if((**cit).status()!=HardBranching::Outgoing) continue;
    // find the partner branchings
    tShowerParticlePtr newPartner=(*cit)->branchingParticle()->partner();
    if(!newPartner) continue;
    tHardBranchingPtr branch;
    for( set<HardBranchingPtr>::iterator clt = branchings.begin();
	 clt != branchings.end(); ++clt ) {
      if(cit==clt) continue;
      if((**clt).branchingParticle()==newPartner) {
	(**cit).colourPartner(*clt);
 	branch=*clt;
	break;
      }
    }
    if((**decay->incoming().begin()).branchingParticle()==newPartner) {
      (**cit).colourPartner(*decay->incoming().begin());
      branch = *decay->incoming().begin();
    }
    // final-state colour partner
    if(branch->status()==HardBranching::Outgoing) {
      Boost boost=((*cit)->pVector()+branch->pVector()).findBoostToCM();
      Lorentz5Momentum pcm = branch->pVector();
      pcm.boost(boost);
      Lorentz5Momentum nvect = Lorentz5Momentum(ZERO,pcm.vect());
      nvect.boost( -boost);
      (*cit)->nVector(nvect);
    }
    // initial-state colour partner
    else {
      Boost boost=branch->pVector().findBoostToCM();
      Lorentz5Momentum pcm = (*cit)->pVector();
      pcm.boost(boost);
      Lorentz5Momentum nvect = Lorentz5Momentum( ZERO, -pcm.vect());
      nvect.boost( -boost);
      (*cit)->nVector(nvect);
    }
  }
  // now compute the new momenta 
  // and calculate the shower variables
  for(cit=branchings.begin();cit!=branchings.end();++cit) {
    if((**cit).status()!=HardBranching::Outgoing) continue;
    LorentzRotation B=LorentzRotation(-boostv);
    LorentzRotation A=LorentzRotation(boostv),R;
    if((*cit)->branchingParticle()==partner) {
      Lorentz5Momentum qnew;
      Energy2 dot=(*cit)->pVector()*(*cit)->nVector();
      double beta = 0.5*((*cit)->branchingParticle()->momentum().m2()
			 -sqr((*cit)->pVector().mass()))/dot;
      qnew=(*cit)->pVector()+beta*(*cit)->nVector();
      qnew.rescaleMass();
      // compute the boost
      R=B*solveBoost(A*qnew,A*(*cit)->branchingParticle()->momentum())*A;
    }
    else {
      Lorentz5Momentum qnew;
      if((*cit)->branchingParticle()->partner()) {
	Energy2 dot=(*cit)->pVector()*(*cit)->nVector();
	double beta = 0.5*((*cit)->branchingParticle()->momentum().m2()
			   -sqr((*cit)->pVector().mass()))/dot;
	qnew=(*cit)->pVector()+beta*(*cit)->nVector();
	qnew.rescaleMass();
      }
      else {
	qnew = (*cit)->pVector();
      }
      // compute the boost
      R=B*solveBoost(A*qnew,A*(*cit)->branchingParticle()->momentum())*A;
    }
    // reconstruct the momenta
    (*cit)->setMomenta(R,1.0,Lorentz5Momentum());
  }
  if(initial) {
    initial->setMomenta(LorentzRotation(),1.0,Lorentz5Momentum());
  }
  return true;
}

double KinematicsReconstructor::
inverseRescalingFactor(vector<Lorentz5Momentum> pout,
		       vector<Energy> mon, Energy roots) const {
  double lambda=1.;
  if(pout.size()==2) { 
    double mu_q1(pout[0].m()/roots), mu_q2(pout[1].m()/roots);
    double mu_p1(mon[0]/roots)     , mu_p2(mon[1]/roots);
    lambda = 
      ((1.+mu_q1+mu_q2)*(1.-mu_q1-mu_q2)*(mu_q1-1.-mu_q2)*(mu_q2-1.-mu_q1))/
      ((1.+mu_p1+mu_p2)*(1.-mu_p1-mu_p2)*(mu_p1-1.-mu_p2)*(mu_p2-1.-mu_p1));
    if(lambda<0.)
      throw Exception() << "Rescaling factor is imaginary in  KinematicsReconstructor::"
			<< "inverseRescalingFactor lambda^2= " << lambda
			<< Exception::eventerror;
    lambda = sqrt(lambda);
  }
  else {
    unsigned int ntry=0;
    // compute magnitudes once for speed
    vector<Energy2> pmag;
    for(unsigned int ix=0;ix<pout.size();++ix) {
      pmag.push_back(pout[ix].vect().mag2());
    }
    // Newton-Raphson for the rescaling
    vector<Energy> root(pout.size());
    do {
      // compute new energies
      Energy sum(ZERO);
      for(unsigned int ix=0;ix<pout.size();++ix) {
	root[ix] = sqrt(pmag[ix]/sqr(lambda)+sqr(mon[ix]));
	sum+=root[ix];
      }
      // if accuracy reached exit
      if(abs(sum/roots-1.)<1e-10) break;
      // use Newton-Raphson to compute new guess for lambda
      Energy numer(ZERO),denom(ZERO);
      for(unsigned int ix=0;ix<pout.size();++ix) {
	numer +=root[ix];
	denom +=pmag[ix]/root[ix];
      }
      numer-=roots;
      double fact = 1.+sqr(lambda)*numer/denom;
      if(fact<0.) fact=0.5;
      lambda *=fact;
      ++ntry;
    }
    while(ntry<100);
  }
  if(std::isnan(lambda))
    throw Exception() << "Rescaling factor is nan in  KinematicsReconstructor::"
		      << "inverseRescalingFactor " 
		      << Exception::eventerror;
  return lambda;
}

bool KinematicsReconstructor::
deconstructGeneralSystem(HardTreePtr tree,
			 ShowerInteraction type) const {
  // extract incoming and outgoing particles
  ColourSingletShower in,out;
  for(set<HardBranchingPtr>::const_iterator it=tree->branchings().begin();
      it!=tree->branchings().end();++it) {
    if((**it).status()==HardBranching::Incoming) in .jets.push_back(*it);
    else                  out.jets.push_back(*it);
  }
  LorentzRotation toRest,fromRest;
  bool applyBoost(false);
  // do the initial-state reconstruction
  deconstructInitialInitialSystem(applyBoost,toRest,fromRest,
				  tree,in.jets,type);
  // do the final-state reconstruction
  deconstructFinalStateSystem(toRest,fromRest,tree,
			      out.jets,type);
  // only at this point that we can be sure all the reference vectors
  // are correct
  for(set<HardBranchingPtr>::const_iterator it=tree->branchings().begin();
      it!=tree->branchings().end();++it) {
    if((**it).status()==HardBranching::Incoming) continue;
    if((**it).branchingParticle()->coloured())
      (**it).setMomenta(LorentzRotation(),1.,Lorentz5Momentum(),false);
  }
  for(set<HardBranchingPtr>::const_iterator it=tree->incoming().begin();
      it!=tree->incoming().end();++it) {
    (**it).setMomenta(LorentzRotation(),1.,Lorentz5Momentum(),false);
  }
  return true;
}

bool KinematicsReconstructor::deconstructHardJets(HardTreePtr tree,
					      ShowerInteraction type) const {
  // inverse of old recon method
  if(_reconopt == 0) {
    return deconstructGeneralSystem(tree,type);
  }
  else if(_reconopt == 1) {
    return deconstructColourSinglets(tree,type);
  }
  else if(_reconopt == 2) {
    throw Exception() << "Inverse reconstruction is not currently supported for ReconstructionOption Colour2 "
		      << "in KinematicsReconstructor::deconstructHardJets(). Please use one of the other options\n"
		      << Exception::runerror;
  }
  else if(_reconopt == 3 || _reconopt == 4 ) {
    return deconstructColourPartner(tree,type);
  }
  else {
    assert(false);
    return false;
  }
}

bool KinematicsReconstructor::
deconstructColourSinglets(HardTreePtr tree,
			  ShowerInteraction type) const {
  // identify the colour singlet systems
  unsigned int nnun(0),nnii(0),nnif(0),nnf(0),nni(0);
  vector<ColourSingletShower> 
    systems(identifySystems(tree->branchings(),nnun,nnii,nnif,nnf,nni));
  // now decide what to do
  LorentzRotation toRest,fromRest;
  bool applyBoost(false);
  bool general(false);
  // initial-initial connection and final-state colour singlet systems
  // Drell-Yan type
  if(nnun==0&&nnii==1&&nnif==0&&nnf>0&&nni==0) {
    // reconstruct initial-initial system
    for(unsigned int ix=0;ix<systems.size();++ix) {
      if(systems[ix].type==II) 
	deconstructInitialInitialSystem(applyBoost,toRest,fromRest,tree,
					systems[ix].jets,type);
    }
    if(type!=ShowerInteraction::QCD) {
      combineFinalState(systems);
      general=false;
    }
  }
  // DIS and VBF type
  else if(nnun==0&&nnii==0&&((nnif==1&&nnf>0&&nni==1)||
			     (nnif==2&&       nni==0))) {
    for(unsigned int ix=0;ix<systems.size();++ix) {
      if(systems[ix].type==IF)
	deconstructInitialFinalSystem(tree,systems[ix].jets,type);
    }
  }
  // e+e- type
  else if(nnun==0&&nnii==0&&nnif==0&&nnf>0&&nni==2) {
    // only FS needed
    // but need to boost to rest frame if QED ISR
    Lorentz5Momentum ptotal;
    for(unsigned int ix=0;ix<systems.size();++ix) {
      if(systems[ix].type==I) 
	ptotal += systems[ix].jets[0]->branchingParticle()->momentum();
    }
    toRest = LorentzRotation(ptotal.findBoostToCM());
    fromRest = toRest;
    fromRest.invert();
    if(type!=ShowerInteraction::QCD) {
      combineFinalState(systems);
      general=false;
    }
  }
  // general type
  else {
    general = true;
  }
  // final-state systems except for general recon
  if(!general) {
    for(unsigned int ix=0;ix<systems.size();++ix) {
      if(systems[ix].type==F)
	deconstructFinalStateSystem(toRest,fromRest,tree,
				    systems[ix].jets,type);
    }
    // only at this point that we can be sure all the reference vectors
    // are correct
    for(set<HardBranchingPtr>::const_iterator it=tree->branchings().begin();
	it!=tree->branchings().end();++it) {
      if((**it).status()==HardBranching::Incoming) continue;
      if((**it).branchingParticle()->coloured())
	(**it).setMomenta(LorentzRotation(),1.,Lorentz5Momentum(),false);
    }
    for(set<HardBranchingPtr>::const_iterator it=tree->incoming().begin();
	it!=tree->incoming().end();++it) {
      (**it).setMomenta(LorentzRotation(),1.,Lorentz5Momentum(),false);
    }
    return true;
  }
  else {
    return deconstructGeneralSystem(tree,type);
  }
  return true;
}

bool KinematicsReconstructor::
deconstructColourPartner(HardTreePtr tree,
			 ShowerInteraction type) const {
  Lorentz5Momentum ptotal;
  HardBranchingPtr emitter;
  ColourSingletShower incomingShower,outgoingShower;
  for(set<HardBranchingPtr>::const_iterator it=tree->branchings().begin();
      it!=tree->branchings().end();++it) { 
    if((**it).status()==HardBranching::Incoming) {
      incomingShower.jets.push_back(*it);
      ptotal += (*it)->branchingParticle()->momentum();
      // check for emitting particle
      if((**it).parent() ) {
	if(!emitter)
	  emitter = *it;
	else
	  throw Exception() << "Only one emitting particle allowed in "
			    << "KinematicsReconstructor::deconstructColourPartner()"
			    << Exception::runerror;
      }
    }
    else if ((**it).status()==HardBranching::Outgoing) {
      outgoingShower.jets.push_back(*it);
      // check for emitting particle
      if(!(**it).children().empty() ) {
	if(!emitter)
	  emitter = *it;
	else
	  throw Exception() << "Only one emitting particle allowed in "
			    << "KinematicsReconstructor::deconstructColourPartner()"
			    << Exception::runerror;
      }
    }
  }
  assert(emitter);
  assert(emitter->colourPartner());
  ColourSingletShower system;
  system.jets.push_back(emitter);
  system.jets.push_back(emitter->colourPartner());
  LorentzRotation toRest,fromRest; 
  bool applyBoost(false);
  // identify the colour singlet system
  if(emitter->status()                  == HardBranching::Outgoing &&
     emitter->colourPartner()->status() == HardBranching::Outgoing ) {
    system.type=F;
    // need to boost to rest frame if QED ISR
    if(  !incomingShower.jets[0]->branchingParticle()->coloured() &&
	 !incomingShower.jets[1]->branchingParticle()->coloured() ) {
      Boost boost = ptotal.findBoostToCM();
      toRest   = LorentzRotation( boost);
      fromRest = LorentzRotation(-boost);
    }
    else
      findInitialBoost(ptotal,ptotal,toRest,fromRest);
    deconstructFinalStateSystem(toRest,fromRest,tree,
				system.jets,type);
  }
  else if (emitter->status()                  == HardBranching::Incoming &&
	   emitter->colourPartner()->status() == HardBranching::Incoming) {
    system.type=II;
    deconstructInitialInitialSystem(applyBoost,toRest,fromRest,tree,system.jets,type);
    // make sure the recoil gets applied
    deconstructFinalStateSystem(toRest,fromRest,tree,
				outgoingShower.jets,type);
  }
  else if ((emitter->status()                  == HardBranching::Outgoing &&
	    emitter->colourPartner()->status() == HardBranching::Incoming ) || 
	   (emitter->status()                  == HardBranching::Incoming &&
	    emitter->colourPartner()->status() == HardBranching::Outgoing)) {
    system.type=IF;
    // enusre incoming first
    if(system.jets[0]->status() == HardBranching::Outgoing)
      swap(system.jets[0],system.jets[1]);
    deconstructInitialFinalSystem(tree,system.jets,type);
  }
  else {
    throw Exception() << "Unknown type of system in "
		      << "KinematicsReconstructor::deconstructColourPartner()"
		      << Exception::runerror;
  }
  // only at this point that we can be sure all the reference vectors
  // are correct
  for(set<HardBranchingPtr>::const_iterator it=tree->branchings().begin();
      it!=tree->branchings().end();++it) {
    if((**it).status()==HardBranching::Incoming) continue;
    if((**it).branchingParticle()->coloured())
      (**it).setMomenta(LorentzRotation(),1.,Lorentz5Momentum(),false);
  }
  for(set<HardBranchingPtr>::const_iterator it=tree->incoming().begin();
      it!=tree->incoming().end();++it) {
    (**it).setMomenta(LorentzRotation(),1.,Lorentz5Momentum(),false);
  }
  for(set<HardBranchingPtr>::const_iterator it=tree->branchings().begin();
      it!=tree->branchings().end();++it) {
    if((**it).status()!=HardBranching::Incoming) continue;
    if(*it==system.jets[0] || *it==system.jets[1]) continue;
    if((**it).branchingParticle()->momentum().z()>ZERO) {
      (**it).z((**it).branchingParticle()->momentum().plus()/(**it).beam()->momentum().plus());
    }
    else {
      (**it).z((**it).branchingParticle()->momentum().minus()/(**it).beam()->momentum().minus());
    }
  }
  return true;
}

void KinematicsReconstructor::
reconstructInitialFinalSystem(vector<ShowerProgenitorPtr> jets) const {
  Lorentz5Momentum pin[2],pout[2],pbeam;
  for(unsigned int ix=0;ix<jets.size();++ix) {
    // final-state parton
    if(jets[ix]->progenitor()->isFinalState()) {
      pout[0] +=jets[ix]->progenitor()->momentum();
      _progenitor = jets[ix]->progenitor();
      if(jets[ix]->reconstructed()==ShowerProgenitor::notReconstructed) {
	reconstructTimeLikeJet(jets[ix]->progenitor());
	jets[ix]->reconstructed(ShowerProgenitor::done);
      }
    }
    // initial-state parton
    else {
      pin[0]  +=jets[ix]->progenitor()->momentum();
      if(jets[ix]->progenitor()->showerKinematics()) {
	pbeam = jets[ix]->progenitor()->showerBasis()->getBasis()[0];
      }
      else {
	if ( jets[ix]->original()->parents().empty() ) {
	  pbeam = jets[ix]->progenitor()->momentum();
	}
	else {
	  pbeam = jets[ix]->original()->parents()[0]->momentum();
	}
      }
      if(jets[ix]->reconstructed()==ShowerProgenitor::notReconstructed) {
	reconstructSpaceLikeJet(jets[ix]->progenitor());
	jets[ix]->reconstructed(ShowerProgenitor::done);
      }
      assert(!jets[ix]->original()->parents().empty());
    }
  }
  // add intrinsic pt if needed
  addIntrinsicPt(jets);
  // momenta after showering
  for(unsigned int ix=0;ix<jets.size();++ix) {
    if(jets[ix]->progenitor()->isFinalState())
      pout[1] += jets[ix]->progenitor()->momentum();
    else
      pin[1]  += jets[ix]->progenitor()->momentum();
  }
  // work out the boost to the Breit frame
  Lorentz5Momentum pa = pout[0]-pin[0];
  Axis axis(pa.vect().unit());
  LorentzRotation rot;
  double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
  if ( sinth > 1.e-9 )
    rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  rot.rotateX(Constants::pi);
  rot.boostZ( pa.e()/pa.vect().mag());
  Lorentz5Momentum ptemp=rot*pbeam;
  Boost trans = -1./ptemp.e()*ptemp.vect();
  trans.setZ(0.);
  if ( trans.mag2() - 1. >= 0. ) throw KinematicsReconstructionVeto();
  rot.boost(trans);
  pa *=rot;
  // project and calculate rescaling
  // reference vectors
  Lorentz5Momentum n1(ZERO,ZERO,-pa.z(),-pa.z());
  Lorentz5Momentum n2(ZERO,ZERO, pa.z(),-pa.z());
  Energy2 n1n2 = n1*n2;
  // decompose the momenta
  Lorentz5Momentum qbp=rot*pin[1],qcp=rot*pout[1];
  qbp.rescaleMass();
  qcp.rescaleMass();
  double a[2],b[2];
  a[0] = n2*qbp/n1n2;
  b[0] = n1*qbp/n1n2;
  Lorentz5Momentum qperp = qbp-a[0]*n1-b[0]*n2;
  b[1] = 0.5;
  a[1] = 0.5*(qcp.m2()-qperp.m2())/n1n2/b[1];
  double kb;
  if(a[0]!=0.) {
    double A(0.5*a[0]),B(b[0]*a[0]-a[1]*b[1]-0.25),C(-0.5*b[0]);
    if(sqr(B)-4.*A*C<0.) throw KinematicsReconstructionVeto();
    kb = 0.5*(-B+sqrt(sqr(B)-4.*A*C))/A;
  }
  else {
    kb = 0.5*b[0]/(b[0]*a[0]-a[1]*b[1]-0.25);
  }
  // changed to improve stability
  if(kb==0.) throw KinematicsReconstructionVeto();
  if ( a[1]>b[1] && abs(a[1]) < 1e-12 )
    throw KinematicsReconstructionVeto();
  if ( a[1]<=b[1] && abs(0.5+b[0]/kb) < 1e-12 )
    throw KinematicsReconstructionVeto();
  double kc = (a[1]>b[1]) ? (a[0]*kb-0.5)/a[1] : b[1]/(0.5+b[0]/kb);
  if(kc==0.) throw KinematicsReconstructionVeto();
  Lorentz5Momentum pnew[2] = { a[0]*kb*n1+b[0]/kb*n2+qperp,
			       a[1]*kc*n1+b[1]/kc*n2+qperp};
  LorentzRotation rotinv=rot.inverse();
  for(unsigned int ix=0;ix<jets.size();++ix) {
    if(jets[ix]->progenitor()->isFinalState()) {
      deepTransform(jets[ix]->progenitor(),rot);
      deepTransform(jets[ix]->progenitor(),solveBoost(pnew[1],qcp));
      Energy delta = jets[ix]->progenitor()->momentum().m()-jets[ix]->progenitor()->momentum().mass();
      if ( abs(delta) > MeV ) throw KinematicsReconstructionVeto();
      deepTransform(jets[ix]->progenitor(),rotinv);
    }
    else {
      tPPtr parent;
      boostChain(jets[ix]->progenitor(),rot,parent);
      boostChain(jets[ix]->progenitor(),solveBoostZ(pnew[0],qbp),parent);
      // check the first boost worked, and if not apply small correction to
      // fix energy/momentum conservation
      // this is a kludge but it reduces momentum non-conservation dramatically
      Lorentz5Momentum pdiff = pnew[0]-jets[ix]->progenitor()->momentum();
      Energy2 delta = sqr(pdiff.x())+sqr(pdiff.y())+sqr(pdiff.z())+sqr(pdiff.t());
      unsigned int ntry=0;
      while(delta>1e-6*GeV2 && ntry<5 ) {
	ntry +=1;
	boostChain(jets[ix]->progenitor(),solveBoostZ(pnew[0],jets[ix]->progenitor()->momentum()),parent);
	pdiff = pnew[0]-jets[ix]->progenitor()->momentum();
	delta = sqr(pdiff.x())+sqr(pdiff.y())+sqr(pdiff.z())+sqr(pdiff.t());
      }
      // apply test in breit-frame
      Lorentz5Momentum ptest1 = parent->momentum();
      Lorentz5Momentum ptest2 = rot*pbeam;
      if(ptest1.z()/ptest2.z()<0. || ptest1.z()/ptest2.z()>1.)
      	throw KinematicsReconstructionVeto();
      boostChain(jets[ix]->progenitor(),rotinv,parent);
    }
  }
}

bool KinematicsReconstructor::addIntrinsicPt(vector<ShowerProgenitorPtr> jets) const {
  bool added=false;
  // add the intrinsic pt if needed
  for(unsigned int ix=0;ix<jets.size();++ix) {
    // only for initial-state particles which haven't radiated
    if(jets[ix]->progenitor()->isFinalState()||
       jets[ix]->hasEmitted()||
       jets[ix]->reconstructed()==ShowerProgenitor::dontReconstruct) continue;
    if(_intrinsic.find(jets[ix])==_intrinsic.end()) continue;
    pair<Energy,double> pt=_intrinsic[jets[ix]];
    Energy etemp = jets[ix]->original()->parents()[0]->momentum().z();
    Lorentz5Momentum 
      p_basis(ZERO, ZERO, etemp, abs(etemp)),
      n_basis(ZERO, ZERO,-etemp, abs(etemp));
    double alpha = jets[ix]->progenitor()->x();
    double beta  = 0.5*(sqr(jets[ix]->progenitor()->data().mass())+
			sqr(pt.first))/alpha/(p_basis*n_basis);
    Lorentz5Momentum pnew=alpha*p_basis+beta*n_basis;
    pnew.setX(pt.first*cos(pt.second));
    pnew.setY(pt.first*sin(pt.second));
    pnew.rescaleMass();
    jets[ix]->progenitor()->set5Momentum(pnew);
    added = true;
  }
  return added;
}

namespace {

double defaultSolveBoostGamma(const double & betam,const Energy2 & kps,
			      const Energy2 & qs, const Energy2 & Q2,
			      const Energy & kp,
			      const Energy & q, const Energy & qE) {
  if(betam<0.5) {
    return 1./sqrt(1.-sqr(betam));
  }
  else {
    return ( kps+ qs + Q2)/
      sqrt(2.*kps*qs + kps*Q2 + qs*Q2 + sqr(Q2) + 2.*q*qE*kp*sqrt(kps + Q2));
  }
}

}

LorentzRotation KinematicsReconstructor::
solveBoost(const double k, const Lorentz5Momentum & newq, 
	   const Lorentz5Momentum & oldp ) const {
  Energy q = newq.vect().mag(); 
  Energy2 qs = sqr(q); 
  Energy2 Q2 = newq.mass2();
  Energy kp = k*(oldp.vect().mag()); 
  Energy2 kps = sqr(kp);
  double betam = (q*newq.e() - kp*sqrt(kps + Q2))/(kps + qs + Q2); 
  if ( abs(betam) - 1. >= 0. ) throw KinematicsReconstructionVeto();
  Boost beta = -betam*(k/kp)*oldp.vect();
  double gamma = 0.;
  if(Q2/sqr(oldp.e())>1e-4) {
    gamma = defaultSolveBoostGamma(betam,kps,qs,Q2,kp,q,newq.e());
  }
  else {
    if(k>0) {
      gamma = 4.*kps*qs/sqr(kps +qs)  + 2.*sqr(kps-qs)*Q2/pow<3,1>(kps +qs) 
  	- 0.25*( sqr(kps) + 14.*kps*qs + sqr(qs))*sqr(kps-qs)/(pow<4,1>(kps +qs)*kps*qs)*sqr(Q2);
    }
    else { 
      gamma = 0.25*sqr(Q2)/(kps*qs)*(1. - 0.5*(kps+qs)/(kps*qs)*Q2);
    }
    if(gamma<=0.) throw KinematicsReconstructionVeto();
    gamma = 1./sqrt(gamma);
    if(gamma>2.) gamma = defaultSolveBoostGamma(betam,kps,qs,Q2,kp,q,newq.e());
  }
  // note that (k/kp)*oldp.vect() = oldp.vect()/oldp.vect().mag() but cheaper.
  ThreeVector<Energy2> ax = newq.vect().cross( oldp.vect() );
  double delta;
  if (newq.x()*oldp.x()+newq.y()*oldp.y()+newq.z()*oldp.z()< 1e-16*GeV2) {
    throw KinematicsReconstructionVeto();
  }else{
    delta = newq.vect().angle( oldp.vect() );
  }
  
  
  LorentzRotation R;
  using Constants::pi;
  Energy2 scale1 = sqr(newq.x())+ sqr(newq.y())+sqr(newq.z());
  Energy2 scale2 = sqr(oldp.x())+ sqr(oldp.y())+sqr(oldp.z());
  if ( ax.mag2()/scale1/scale2 > 1e-28 ) {
    R.rotate( delta, unitVector(ax) ).boost( beta , gamma );
  } 
  else if(abs(delta-pi)/pi < 0.001) {
    double phi=2.*pi*UseRandom::rnd();
    Axis axis(cos(phi),sin(phi),0.);
    axis.rotateUz(newq.vect().unit());
    R.rotate(delta,axis).boost( beta , gamma );
  }
  else {
    R.boost( beta , gamma );
  }
  return R;
}

LorentzRotation KinematicsReconstructor::solveBoost(const Lorentz5Momentum & q, 
						const Lorentz5Momentum & p ) const {
  Energy modp = p.vect().mag();
  Energy modq = q.vect().mag();
  double betam = (p.e()*modp-q.e()*modq)/(sqr(modq)+sqr(modp)+p.mass2());
  if ( abs(betam)-1. >= 0. ) throw KinematicsReconstructionVeto();
  Boost beta = -betam*q.vect().unit();
  ThreeVector<Energy2> ax = p.vect().cross( q.vect() ); 
  double delta = p.vect().angle( q.vect() );
  LorentzRotation R;
  using Constants::pi;
  if ( beta.mag2() - 1. >= 0. ) throw KinematicsReconstructionVeto();
  if ( ax.mag2()/GeV2/MeV2 > 1e-16 ) {
    R.rotate( delta, unitVector(ax) ).boost( beta );
  } 
  else {
    R.boost( beta );
  } 
  return R;
}

LorentzRotation KinematicsReconstructor::solveBoostZ(const Lorentz5Momentum & q, 
						 const Lorentz5Momentum & p ) const {
  static const double eps = 1e-6;
  LorentzRotation R;
  double beta;
  Energy2 mt2  = p.mass()<ZERO ? -sqr(p.mass())+sqr(p.x())+sqr(p.y()) : sqr(p.mass())+sqr(p.x())+sqr(p.y())  ;
  double ratio = mt2/(sqr(p.t())+sqr(q.t()));
  if(abs(ratio)>eps) {
    double erat  = (q.t()+q.z())/(p.t()+p.z());
    Energy2 den = mt2*(erat+1./erat);
    Energy2 num = (q.z()-p.z())*(q.t()+p.t()) + (p.z()+q.z())*(p.t()-q.t());
    beta = num/den;
    if ( abs(beta) - 1. >= 0. ) throw KinematicsReconstructionVeto();
    R.boostZ(beta);
  }
  else {
    double er = sqr(p.t()/q.t());
    double x = ratio+0.125*(er+10.+1./er)*sqr(ratio);
    beta = -(p.t()-q.t())*(p.t()+q.t())/(sqr(p.t())+sqr(q.t()))*(1.+x);
    double gamma = (4.*sqr(p.t()*q.t()) +sqr(p.t()-q.t())*sqr(p.t()+q.t())*
		    (-2.*x+sqr(x)))/sqr(sqr(p.t())+sqr(q.t()));
    if ( abs(beta) - 1. >= 0. ) throw KinematicsReconstructionVeto();
    gamma  = 1./sqrt(gamma);
    R.boost(0.,0.,beta,gamma);
  }
  Lorentz5Momentum ptest = R*p;
  if(ptest.z()/q.z() < 0. || ptest.t()/q.t() < 0. ) {
    throw KinematicsReconstructionVeto();
  }
  return R;
}

void KinematicsReconstructor::
reconstructFinalStateSystem(bool applyBoost, 
			    const LorentzRotation &   toRest,
			    const LorentzRotation & fromRest, 
			    vector<ShowerProgenitorPtr> jets) const {
  LorentzRotation trans = applyBoost? toRest : LorentzRotation();
  // special for case of individual particle
  if(jets.size()==1) {
    deepTransform(jets[0]->progenitor(),trans);
    deepTransform(jets[0]->progenitor(),fromRest);
    return;
  }
  bool radiated(false);
  // find the hard process centre-of-mass energy
  Lorentz5Momentum pcm;
  // check if radiated and calculate total momentum
  for(unsigned int ix=0;ix<jets.size();++ix) {
    radiated |=jets[ix]->hasEmitted();
    pcm += jets[ix]->progenitor()->momentum();
  }
  if(applyBoost) pcm *= trans;
  // check if in CMF frame
  Boost beta_cm = pcm.findBoostToCM();
  bool gottaBoost(false);
  if(beta_cm.mag() > 1e-12) {
    gottaBoost = true;
    trans.boost(beta_cm);
  }
  // collection of pointers to initial hard particle and jet momenta
  // for final boosts
  JetKinVect jetKinematics;
  vector<ShowerProgenitorPtr>::const_iterator cit;
  for(cit = jets.begin(); cit != jets.end(); cit++) {
    JetKinStruct tempJetKin;      
    tempJetKin.parent = (*cit)->progenitor();
    if(applyBoost || gottaBoost) {
      deepTransform(tempJetKin.parent,trans);
    }
    tempJetKin.p = (*cit)->progenitor()->momentum();
    _progenitor=tempJetKin.parent;
    if((**cit).reconstructed()==ShowerProgenitor::notReconstructed) {
      radiated |= reconstructTimeLikeJet((*cit)->progenitor());
      (**cit).reconstructed(ShowerProgenitor::done);
    }
    else {
      radiated |= !(*cit)->progenitor()->children().empty();
    }
    tempJetKin.q = (*cit)->progenitor()->momentum();
    jetKinematics.push_back(tempJetKin);
  }
  if(_finalFinalWeight && jetKinematics.size()==2) {
    Energy m1 = jetKinematics[0].q.m();
    Energy m2 = jetKinematics[1].q.m();
    Energy m0 = pcm.m();
    if(m0<m1+m2)  throw KinematicsReconstructionVeto();
    Energy4 lambdaNew = (sqr(m0)-sqr(m1-m2))*(sqr(m0)-sqr(m1+m2));
    m1 = jetKinematics[0].p.m();
    m2 = jetKinematics[1].p.m();
    Energy4 lambdaOld = (sqr(m0)-sqr(m1-m2))*(sqr(m0)-sqr(m1+m2));
    if(UseRandom::rnd()>sqrt(lambdaNew/lambdaOld))
      throw KinematicsReconstructionVeto();
  }
  // default option rescale everything with the same factor
  // find the rescaling factor
  double k = 0.0;
  if(radiated) {
    k = solveKfactor(pcm.m(), jetKinematics);
    // perform the rescaling and boosts
    for(JetKinVect::iterator it = jetKinematics.begin();
	it != jetKinematics.end(); ++it) {
      LorentzRotation Trafo = solveBoost(k, it->q, it->p);
      deepTransform(it->parent,Trafo);
    }
  }
  // apply the final boosts
  if(gottaBoost || applyBoost) { 
    LorentzRotation finalBoosts;
    if(gottaBoost) finalBoosts.boost(-beta_cm);
    if(applyBoost) finalBoosts.transform(fromRest);
    for(JetKinVect::iterator it = jetKinematics.begin();
	it != jetKinematics.end(); ++it) {
      deepTransform(it->parent,finalBoosts);
    }
  }
}

void KinematicsReconstructor::
reconstructInitialInitialSystem(bool & applyBoost, 
				LorentzRotation &   toRest, 
				LorentzRotation & fromRest,  
				vector<ShowerProgenitorPtr> jets) const {
  bool radiated = false;
  Lorentz5Momentum pcm;
  // check whether particles radiated and calculate total momentum
  for( unsigned int ix = 0; ix < jets.size(); ++ix ) {
    radiated |= jets[ix]->hasEmitted();
    pcm += jets[ix]->progenitor()->momentum();
    if(jets[ix]->original()->parents().empty()) return;
  }
  pcm.rescaleMass();
  // check if intrinsic pt to be added
  radiated |= !_intrinsic.empty();
  // if no radiation return
  if(!radiated) {
    for(unsigned int ix=0;ix<jets.size();++ix) {
      if(jets[ix]->reconstructed()==ShowerProgenitor::notReconstructed)
	jets[ix]->reconstructed(ShowerProgenitor::done);
    }
    return;
  }
  // initial state shuffling
  applyBoost=false;
  vector<Lorentz5Momentum> p, pq, p_in;
  vector<Energy> pts;
  for(unsigned int ix=0;ix<jets.size();++ix) {
    // add momentum to vector
    p_in.push_back(jets[ix]->progenitor()->momentum());
    // reconstruct the jet
    if(jets[ix]->reconstructed()==ShowerProgenitor::notReconstructed) {
      radiated |= reconstructSpaceLikeJet(jets[ix]->progenitor());
      jets[ix]->reconstructed(ShowerProgenitor::done);
    }
    assert(!jets[ix]->original()->parents().empty());
    Energy etemp = jets[ix]->original()->parents()[0]->momentum().z();
    Lorentz5Momentum ptemp = Lorentz5Momentum(ZERO, ZERO, etemp, abs(etemp));
    pq.push_back(ptemp);
    pts.push_back(jets[ix]->highestpT());
  }
  // add the intrinsic pt if needed
  radiated |=addIntrinsicPt(jets);
  for(unsigned int ix=0;ix<jets.size();++ix) {
    p.push_back(jets[ix]->progenitor()->momentum());
  }
  double x1 = p_in[0].z()/pq[0].z();
  double x2 = p_in[1].z()/pq[1].z();
  vector<double> beta=initialStateRescaling(x1,x2,p_in[0]+p_in[1],p,pq,pts);
  // if not need don't apply boosts
  if(!(radiated && p.size() == 2 && pq.size() == 2)) return;
  applyBoost=true;
  // apply the boosts
  Lorentz5Momentum newcmf;
  for(unsigned int ix=0;ix<jets.size();++ix) {
    tPPtr toBoost = jets[ix]->progenitor();
    Boost betaboost(0, 0, beta[ix]);
    tPPtr parent;
    boostChain(toBoost, LorentzRotation(0.,0.,beta[ix]),parent);
    if(parent->momentum().e()/pq[ix].e()>1.||
       parent->momentum().z()/pq[ix].z()>1.) throw KinematicsReconstructionVeto();
    newcmf+=toBoost->momentum();
  }
  if(newcmf.m()<ZERO||newcmf.e()<ZERO) throw KinematicsReconstructionVeto();
  findInitialBoost(pcm,newcmf,toRest,fromRest);
}

void KinematicsReconstructor::
deconstructInitialInitialSystem(bool & applyBoost,
				LorentzRotation & toRest,
				LorentzRotation & fromRest,
				HardTreePtr tree,
				vector<HardBranchingPtr> jets,
				ShowerInteraction) const {
  assert(jets.size()==2);
  // put beam with +z first
  if(jets[0]->beam()->momentum().z()<ZERO) swap(jets[0],jets[1]);
  // get the momenta of the particles
  vector<Lorentz5Momentum> pin,pq;
  for(unsigned int ix=0;ix<jets.size();++ix) {
    pin.push_back(jets[ix]->branchingParticle()->momentum());
    Energy etemp = jets[ix]->beam()->momentum().z();
    pq.push_back(Lorentz5Momentum(ZERO, ZERO,etemp, abs(etemp)));
  }
  // calculate the rescaling
  double x[2];
  Lorentz5Momentum pcm=pin[0]+pin[1];
  assert(pcm.mass2()>ZERO);
  pcm.rescaleMass();
  vector<double> boost = inverseInitialStateRescaling(x[0],x[1],pcm,pin,pq);
  set<HardBranchingPtr>::const_iterator cjt=tree->incoming().begin();
  HardBranchingPtr incoming[2];
  incoming[0] = *cjt;
  ++cjt;
  incoming[1] = *cjt;
  if((*tree->incoming().begin())->beam()->momentum().z()/pq[0].z()<0.)
    swap(incoming[0],incoming[1]);
  // apply the boost the the particles
  unsigned int iswap[2]={1,0};
  for(unsigned int ix=0;ix<2;++ix) {
    LorentzRotation R(0.,0.,-boost[ix]);
    incoming[ix]->pVector(pq[ix]);
    incoming[ix]->nVector(pq[iswap[ix]]);
    incoming[ix]->setMomenta(R,1.,Lorentz5Momentum());
    jets[ix]->showerMomentum(x[ix]*jets[ix]->pVector());
  }
  // and calculate the boosts 
  applyBoost=true;
  // do one boost
  if(_initialBoost==0) {
    toRest   = LorentzRotation(-pcm.boostVector());
  }
  else if(_initialBoost==1) {
    // first the transverse boost
    Energy pT = sqrt(sqr(pcm.x())+sqr(pcm.y()));
    double beta = -pT/pcm.t();
    toRest=LorentzRotation(Boost(beta*pcm.x()/pT,beta*pcm.y()/pT,0.));
    // the longitudinal 
    beta = pcm.z()/sqrt(pcm.m2()+sqr(pcm.z()));
    toRest.boost(Boost(0.,0.,-beta));
  }
  else
    assert(false);
  fromRest = LorentzRotation((jets[0]->showerMomentum()+
  			      jets[1]->showerMomentum()).boostVector());
}

void KinematicsReconstructor::
deconstructFinalStateSystem(const LorentzRotation &   toRest,
			    const LorentzRotation & fromRest,
			    HardTreePtr tree, vector<HardBranchingPtr> jets,
			    ShowerInteraction type) const {
  LorentzRotation trans = toRest;
  if(jets.size()==1) {
    Lorentz5Momentum pnew = toRest*(jets[0]->branchingParticle()->momentum());
    pnew *= fromRest;
    jets[0]->      original(pnew);
    jets[0]->showerMomentum(pnew);
    // find the colour partners
    ShowerParticleVector particles;
    vector<Lorentz5Momentum> ptemp;
    set<HardBranchingPtr>::const_iterator cjt;
    for(cjt=tree->branchings().begin();cjt!=tree->branchings().end();++cjt) {
      ptemp.push_back((**cjt).branchingParticle()->momentum());
      (**cjt).branchingParticle()->set5Momentum((**cjt).showerMomentum());
      particles.push_back((**cjt).branchingParticle());
    }
    dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->partnerFinder()
      ->setInitialEvolutionScales(particles,false,type,false);
    // calculate the reference vectors
    unsigned int iloc(0);
    set<HardBranchingPtr>::iterator clt;
    for(cjt=tree->branchings().begin();cjt!=tree->branchings().end();++cjt) {
      // reset the momentum
      (**cjt).branchingParticle()->set5Momentum(ptemp[iloc]);
      ++iloc;
      // sort out the partners
      tShowerParticlePtr partner = 
	(*cjt)->branchingParticle()->partner();
      if(!partner) continue;
      for(clt=tree->branchings().begin();clt!=tree->branchings().end();++clt) {
	if((**clt).branchingParticle()==partner) {
	  (**cjt).colourPartner(*clt);
	  break;
	}
      }
      tHardBranchingPtr branch;
      for(clt=tree->branchings().begin();clt!=tree->branchings().end();++clt) {
	if(clt==cjt) continue;
	if((*clt)->branchingParticle()==partner) {
	  branch=*clt;
	  break;
	}
      }
    }
    return;
  }
  vector<HardBranchingPtr>::iterator cit;
  vector<Lorentz5Momentum> pout;
  vector<Energy> mon;
  Lorentz5Momentum pin;
  for(cit=jets.begin();cit!=jets.end();++cit) {
    pout.push_back((*cit)->branchingParticle()->momentum());
    mon.push_back(findMass(*cit));
    pin+=pout.back();
  }
  // boost all the momenta to the rest frame of the decaying particle
  pin.rescaleMass();
  pin *=trans;
  Boost beta_cm = pin.findBoostToCM();
  bool gottaBoost(false);
  if(beta_cm.mag() > 1e-12) {
    gottaBoost = true;
    trans.boost(beta_cm);
    pin.boost(beta_cm);
  }
  for(unsigned int ix=0;ix<pout.size();++ix) {
    pout[ix].transform(trans);
  }
  // rescaling factor
  double lambda=inverseRescalingFactor(pout,mon,pin.mass());
  if (lambda< 1.e-10) throw KinematicsReconstructionVeto();
  // now calculate the p reference vectors 
  for(unsigned int ix=0;ix<jets.size();++ix) {
    Lorentz5Momentum pvect = jets[ix]->branchingParticle()->momentum();
    pvect.transform(trans);
    pvect /= lambda;
    pvect.setMass(mon[ix]);
    pvect.rescaleEnergy();
    if(gottaBoost) pvect.boost(-beta_cm);
    pvect.transform(fromRest);
    jets[ix]->pVector(pvect);
    jets[ix]->showerMomentum(pvect);
  }
  // find the colour partners
  ShowerParticleVector particles;
  vector<Lorentz5Momentum> ptemp;
  set<HardBranchingPtr>::const_iterator cjt;
  for(cjt=tree->branchings().begin();cjt!=tree->branchings().end();++cjt) {
    ptemp.push_back((**cjt).branchingParticle()->momentum());
    (**cjt).branchingParticle()->set5Momentum((**cjt).showerMomentum());
    particles.push_back((**cjt).branchingParticle());
  }
  dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->partnerFinder()
    ->setInitialEvolutionScales(particles,false,type,false);
  // calculate the reference vectors
  unsigned int iloc(0);
  set<HardBranchingPtr>::iterator clt;
  for(cjt=tree->branchings().begin();cjt!=tree->branchings().end();++cjt) {
    // reset the momentum
    (**cjt).branchingParticle()->set5Momentum(ptemp[iloc]);
    ++iloc;
  }
  for(cjt=tree->branchings().begin();cjt!=tree->branchings().end();++cjt) {
    // sort out the partners
    tShowerParticlePtr partner = 
      (*cjt)->branchingParticle()->partner();
    if(!partner) continue;
    for(clt=tree->branchings().begin();clt!=tree->branchings().end();++clt) {
      if((**clt).branchingParticle()==partner) {
	(**cjt).colourPartner(*clt);
	break;
      }
    }
    tHardBranchingPtr branch;
    for(clt=tree->branchings().begin();clt!=tree->branchings().end();++clt) {
      if(clt==cjt) continue;
      if((*clt)->branchingParticle()==partner) {
 	branch=*clt;
 	break;
      }
    }
    // compute the reference vectors
    // both incoming, should all ready be done
    if((**cjt).status()==HardBranching::Incoming &&
       (**clt).status()==HardBranching::Incoming) {
      continue;
    }
    // both outgoing
    else if((**cjt).status()!=HardBranching::Incoming&&
	    branch->status()==HardBranching::Outgoing) {
      Boost boost=((*cjt)->pVector()+branch->pVector()).findBoostToCM();
      Lorentz5Momentum pcm = branch->pVector();
      pcm.boost(boost);
      Lorentz5Momentum nvect = Lorentz5Momentum(ZERO,pcm.vect());
      nvect.boost( -boost);
      (**cjt).nVector(nvect);
    }
    else if((**cjt).status()==HardBranching::Incoming) {
      Lorentz5Momentum pa = -(**cjt).showerMomentum()+branch->showerMomentum();
      Lorentz5Momentum pb =  (**cjt).showerMomentum();
      Axis axis(pa.vect().unit());
      LorentzRotation rot;
      double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
      rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      rot.rotateX(Constants::pi);
      rot.boostZ( pa.e()/pa.vect().mag());
      pb*=rot;
      Boost trans = -1./pb.e()*pb.vect();
      trans.setZ(0.);
      rot.boost(trans);
      Energy scale=(**cjt).beam()->momentum().e();
      Lorentz5Momentum pbasis(ZERO,(**cjt).beam()->momentum().vect().unit()*scale);
      Lorentz5Momentum pcm = rot*pbasis;
      rot.invert();
      (**cjt).nVector(rot*Lorentz5Momentum(ZERO,-pcm.vect()));
      tHardBranchingPtr branch2 = *cjt;;      
      while (branch2->parent()) {
	branch2=branch2->parent();
	branch2->nVector(rot*Lorentz5Momentum(ZERO,-pcm.vect()));
      }
    }
    else if(branch->status()==HardBranching::Incoming) {
      (**cjt).nVector(Lorentz5Momentum(ZERO,branch->showerMomentum().vect()));
    }
  }
  // now compute the new momenta 
  for(cjt=tree->branchings().begin();cjt!=tree->branchings().end();++cjt) {
    if(!(*cjt)->branchingParticle()->isFinalState()) continue;
    Lorentz5Momentum qnew;
    if((*cjt)->branchingParticle()->partner()) {
      Energy2 dot=(*cjt)->pVector()*(*cjt)->nVector();
      double beta = 0.5*((*cjt)->branchingParticle()->momentum().m2()
			 -sqr((*cjt)->pVector().mass()))/dot;
      qnew=(*cjt)->pVector()+beta*(*cjt)->nVector();
      qnew.rescaleMass();
    }
    else {
      qnew = (*cjt)->pVector();
    }
    // qnew is the unshuffled momentum in the rest frame of the p basis vectors,
    // for the simple case Z->q qbar g this was checked against analytic formulae.
    // compute the boost
    LorentzRotation R=solveBoost(qnew,
				 toRest*(*cjt)->branchingParticle()->momentum())*toRest;
    (*cjt)->setMomenta(R,1.0,Lorentz5Momentum());  
  }
}

Energy KinematicsReconstructor::momConsEq(double k, 
				      const Energy & root_s, 
				      const JetKinVect & jets) const {
  static const Energy2 eps=1e-8*GeV2;
  Energy dum = ZERO;
  for(JetKinVect::const_iterator it = jets.begin(); it != jets.end(); ++it) {
    Energy2 dum2 = (it->q).m2() + sqr(k)*(it->p).vect().mag2();
    if(dum2 < ZERO) {
      if(dum2 < -eps) throw KinematicsReconstructionVeto();
      dum2 = ZERO;
    }
    dum += sqrt(dum2);
  }
  return dum - root_s; 
}

void KinematicsReconstructor::boostChain(tPPtr p, const LorentzRotation &bv,
				     tPPtr & parent) const {
  if(!p->parents().empty()) boostChain(p->parents()[0], bv,parent);
  else parent=p;
  p->transform(bv);
  if(p->children().size()==2) {
    if(dynamic_ptr_cast<ShowerParticlePtr>(p->children()[1]))
      deepTransform(p->children()[1],bv);
  }
}

namespace {

bool sortJets(ShowerProgenitorPtr j1, ShowerProgenitorPtr j2) {
  return j1->highestpT()>j2->highestpT();
}

}

void KinematicsReconstructor::
reconstructGeneralSystem(vector<ShowerProgenitorPtr> & ShowerHardJets) const {
  // find initial- and final-state systems
  ColourSingletSystem in,out;
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    if(ShowerHardJets[ix]->progenitor()->isFinalState())
      out.jets.push_back(ShowerHardJets[ix]);
    else
      in.jets.push_back(ShowerHardJets[ix]);
  }
  // reconstruct initial-initial system
  LorentzRotation toRest,fromRest;
  bool applyBoost(false);
  // reconstruct initial-initial system
  reconstructInitialInitialSystem(applyBoost,toRest,fromRest,in.jets);
  // reconstruct the final-state systems
  reconstructFinalStateSystem(applyBoost,toRest,fromRest,out.jets);
}


void KinematicsReconstructor::
reconstructFinalFirst(vector<ShowerProgenitorPtr> & ShowerHardJets) const {
  static const Energy2 minQ2 = 1e-4*GeV2;
  map<ShowerProgenitorPtr,bool> used;
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    used[ShowerHardJets[ix]] = false;
  }  // first to the final-state reconstruction of any systems which need it
  set<ShowerProgenitorPtr> outgoing;
  // first find any particles with final state partners
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    if(ShowerHardJets[ix]->progenitor()->isFinalState()&&
       ShowerHardJets[ix]->progenitor()->partner()&&
       ShowerHardJets[ix]->progenitor()->partner()->isFinalState()) outgoing.insert(ShowerHardJets[ix]);
  }
  // then find the colour partners
  if(!outgoing.empty()) {
    set<ShowerProgenitorPtr> partners;
    for(set<ShowerProgenitorPtr>::const_iterator it=outgoing.begin();it!=outgoing.end();++it) {
      for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
	if((**it).progenitor()->partner()==ShowerHardJets[ix]->progenitor()) {
	  partners.insert(ShowerHardJets[ix]);
	  break;
	}
      }
    }
    outgoing.insert(partners.begin(),partners.end());
  }
  // do the final-state reconstruction if needed
  if(!outgoing.empty()) {
    assert(outgoing.size()!=1);
    LorentzRotation toRest,fromRest;
    vector<ShowerProgenitorPtr> outgoingJets(outgoing.begin(),outgoing.end());
    reconstructFinalStateSystem(false,toRest,fromRest,outgoingJets);
  }
  // Now do any initial-final systems which are needed
  vector<ColourSingletSystem> IFSystems;
  // find the systems N.B. can have duplicates
  // find initial-state with FS partners or FS with IS partners
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    if(!ShowerHardJets[ix]->progenitor()->isFinalState()&&
       ShowerHardJets[ix]->progenitor()->partner()&&
       ShowerHardJets[ix]->progenitor()->partner()->isFinalState()) {
      IFSystems.push_back(ColourSingletSystem(IF,ShowerHardJets[ix]));
    }
    else if(ShowerHardJets[ix]->progenitor()->isFinalState()&&
	    ShowerHardJets[ix]->progenitor()->partner()&&
	    !ShowerHardJets[ix]->progenitor()->partner()->isFinalState()) {
      IFSystems.push_back(ColourSingletSystem(IF,ShowerHardJets[ix]));
    }
  }
  // then add the partners
  for(unsigned int is=0;is<IFSystems.size();++is) {
    for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
      if(IFSystems[is].jets[0]->progenitor()->partner()==ShowerHardJets[ix]->progenitor()) {
	IFSystems[is].jets.push_back(ShowerHardJets[ix]);
      }
    }
    // ensure incoming first
    if(IFSystems[is].jets[0]->progenitor()->isFinalState())
      swap(IFSystems[is].jets[0],IFSystems[is].jets[1]);
  }
  if(!IFSystems.empty()) {
    unsigned int istart = UseRandom::irnd(IFSystems.size());
    unsigned int istop=IFSystems.size();
    for(unsigned int is=istart;is<=istop;++is) {
      if(is==IFSystems.size()) {
	if(istart!=0) {
	  istop = istart-1;
	  is=0;
	}
	else break;
      }
      // skip duplicates
      if(used[IFSystems[is].jets[0]] &&
	 used[IFSystems[is].jets[1]] ) continue;
      if(IFSystems[is].jets[0]->original()&&IFSystems[is].jets[0]->original()->parents().empty()) continue;
      Lorentz5Momentum psum;
      for(unsigned int ix=0;ix<IFSystems[is].jets.size();++ix) {
	if(IFSystems[is].jets[ix]->progenitor()->isFinalState())
	  psum += IFSystems[is].jets[ix]->progenitor()->momentum();
	else
	  psum -= IFSystems[is].jets[ix]->progenitor()->momentum();
      }
      if(-psum.m2()>minQ2) {
	reconstructInitialFinalSystem(IFSystems[is].jets);
	for(unsigned int ix=0;ix<IFSystems[is].jets.size();++ix) {
	  used[IFSystems[is].jets[ix]] = true;
	}
      }
    }
  }
  // now we finally need to handle the initial state system
  ColourSingletSystem in,out;
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    if(ShowerHardJets[ix]->progenitor()->isFinalState())
      out.jets.push_back(ShowerHardJets[ix]);
    else
      in.jets.push_back(ShowerHardJets[ix]);
  }
  // reconstruct initial-initial system
  bool doRecon = false;
  for(unsigned int ix=0;ix<in.jets.size();++ix) {
    if(!used[in.jets[ix]]) {
      doRecon = true;
      break;
    }
  }
  LorentzRotation toRest,fromRest;
  bool applyBoost(false);
  if(doRecon) {
    reconstructInitialInitialSystem(applyBoost,toRest,fromRest,in.jets);
  }
  // reconstruct the final-state systems
  if(!doRecon) {
    for(unsigned int ix=0;ix<out.jets.size();++ix) {
      if(!used[out.jets[ix]]) {
	doRecon = true;
	break;
      }
    }
  }
  if(doRecon) {
    reconstructFinalStateSystem(applyBoost,toRest,fromRest,out.jets);
  }
}

void KinematicsReconstructor::
reconstructColourPartner(vector<ShowerProgenitorPtr> & ShowerHardJets) const {
  static const Energy2 minQ2 = 1e-4*GeV2;
  // sort the vector by hardness of emission
  std::sort(ShowerHardJets.begin(),ShowerHardJets.end(),sortJets);
  // map between particles and progenitors for easy lookup
  map<ShowerParticlePtr,ShowerProgenitorPtr> progenitorMap;
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    progenitorMap[ShowerHardJets[ix]->progenitor()] = ShowerHardJets[ix];
  }
  // check that the IF systems can be reconstructed
  bool canReconstruct = true;
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    tShowerParticlePtr progenitor = ShowerHardJets[ix]->progenitor();
    tShowerParticlePtr partner    = progenitor->partner();
    if(!partner) continue;
    else if((progenitor->isFinalState() &&
	     !partner->isFinalState()) ||
	    (!progenitor->isFinalState() &&
	     partner->isFinalState()) ) {
      vector<ShowerProgenitorPtr> jets(2);
      jets[0] = ShowerHardJets[ix];
      jets[1] = progenitorMap[partner];
      Lorentz5Momentum psum;
      for(unsigned int iy=0;iy<jets.size();++iy) {
	if(jets[iy]->progenitor()->isFinalState())
	  psum += jets[iy]->progenitor()->momentum();
	else
	  psum -= jets[iy]->progenitor()->momentum();
      }
      if(-psum.m2()<minQ2) {
	canReconstruct  = false;
	break;
      }
    }
  }
  if(!canReconstruct) {
    reconstructGeneralSystem(ShowerHardJets);
    return;
  }
  map<ShowerProgenitorPtr,bool> used;
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    used[ShowerHardJets[ix]] = false;
  }
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    // skip jets which have already been handled
    if(ShowerHardJets[ix]->reconstructed()==ShowerProgenitor::done) continue;
    // already reconstructed
    if(used[ShowerHardJets[ix]]) continue; 
    // no partner continue
    tShowerParticlePtr progenitor = ShowerHardJets[ix]->progenitor();
    tShowerParticlePtr partner    = progenitor->partner();
    if(!partner) {
      // check if there's a daughter tree which also needs boosting
      Lorentz5Momentum porig = progenitor->momentum();
      map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator tit;
      for(tit  = _currentTree->treelinks().begin();
	  tit != _currentTree->treelinks().end();++tit) {
	// if there is, boost it
	if(tit->second.first && tit->second.second==progenitor) {
	  Lorentz5Momentum pnew = tit->first->incomingLines().begin()
	    ->first->progenitor()->momentum();
	  pnew *=  tit->first->transform();
	  Lorentz5Momentum pdiff = porig-pnew;
	  Energy2 test = sqr(pdiff.x()) + sqr(pdiff.y()) + 
	    sqr(pdiff.z()) + sqr(pdiff.t());
	  LorentzRotation rot;
	  if(test>1e-6*GeV2) rot = solveBoost(porig,pnew);
	  tit->first->transform(rot,false);
	  _treeBoosts[tit->first].push_back(rot);
	}
      }
      ShowerHardJets[ix]->reconstructed(ShowerProgenitor::done);
      continue;
    }
    // do the reconstruction
    // final-final
    if(progenitor->isFinalState() &&
       partner->isFinalState() ) {
      LorentzRotation toRest,fromRest;
      vector<ShowerProgenitorPtr> jets(2);
      jets[0] = ShowerHardJets[ix];
      jets[1] = progenitorMap[partner];
      if(_reconopt==4 && jets[1]->reconstructed()==ShowerProgenitor::notReconstructed)
	jets[1]->reconstructed(ShowerProgenitor::dontReconstruct);
      reconstructFinalStateSystem(false,toRest,fromRest,jets);
      if(_reconopt==4 && jets[1]->reconstructed()==ShowerProgenitor::dontReconstruct)
	jets[1]->reconstructed(ShowerProgenitor::notReconstructed);
      used[jets[0]] = true;
      if(_reconopt==3) used[jets[1]] = true;
    }
    // initial-final
    else if((progenitor->isFinalState() &&
	     !partner->isFinalState()) ||
	    (!progenitor->isFinalState() &&
	     partner->isFinalState()) ) {
      vector<ShowerProgenitorPtr> jets(2);
      jets[0] = ShowerHardJets[ix];
      jets[1] = progenitorMap[partner];
      if(jets[0]->progenitor()->isFinalState()) swap(jets[0],jets[1]);
      if(jets[0]->original()&&jets[0]->original()->parents().empty()) continue;
      Lorentz5Momentum psum;
      for(unsigned int iy=0;iy<jets.size();++iy) {
	if(jets[iy]->progenitor()->isFinalState())
	  psum += jets[iy]->progenitor()->momentum();
	else
	  psum -= jets[iy]->progenitor()->momentum();
      }
      if(_reconopt==4 && progenitorMap[partner]->reconstructed()==ShowerProgenitor::notReconstructed)
	progenitorMap[partner]->reconstructed(ShowerProgenitor::dontReconstruct);
      reconstructInitialFinalSystem(jets);
      if(_reconopt==4 && progenitorMap[partner]->reconstructed()==ShowerProgenitor::dontReconstruct)
	progenitorMap[partner]->reconstructed(ShowerProgenitor::notReconstructed);
      used[ShowerHardJets[ix]] = true;
      if(_reconopt==3) used[progenitorMap[partner]] = true;
    }
    // initial-initial
    else if(!progenitor->isFinalState() &&
	    !partner->isFinalState() ) {
      ColourSingletSystem in,out;
      in.jets.push_back(ShowerHardJets[ix]);
      in.jets.push_back(progenitorMap[partner]);
      for(unsigned int iy=0;iy<ShowerHardJets.size();++iy) {
	if(ShowerHardJets[iy]->progenitor()->isFinalState())
	  out.jets.push_back(ShowerHardJets[iy]);
      }
      LorentzRotation toRest,fromRest;
      bool applyBoost(false);
      if(_reconopt==4 && in.jets[1]->reconstructed()==ShowerProgenitor::notReconstructed)
	in.jets[1]->reconstructed(ShowerProgenitor::dontReconstruct);
      reconstructInitialInitialSystem(applyBoost,toRest,fromRest,in.jets);
      if(_reconopt==4 && in.jets[1]->reconstructed()==ShowerProgenitor::dontReconstruct)
	in.jets[1]->reconstructed(ShowerProgenitor::notReconstructed);
      used[in.jets[0]] = true;
      if(_reconopt==3) used[in.jets[1]] = true;
      for(unsigned int iy=0;iy<out.jets.size();++iy) {
	if(out.jets[iy]->reconstructed()==ShowerProgenitor::notReconstructed)
	  out.jets[iy]->reconstructed(ShowerProgenitor::dontReconstruct);
      }
      // reconstruct the final-state systems
      LorentzRotation finalBoosts;
      finalBoosts.transform(  toRest);
      finalBoosts.transform(fromRest);
      for(unsigned int iy=0;iy<out.jets.size();++iy) {
	deepTransform(out.jets[iy]->progenitor(),finalBoosts);
      }
      for(unsigned int iy=0;iy<out.jets.size();++iy) {
	if(out.jets[iy]->reconstructed()==ShowerProgenitor::dontReconstruct)
	  out.jets[iy]->reconstructed(ShowerProgenitor::notReconstructed);
      }
    }
  }
}

bool KinematicsReconstructor::
inverseDecayRescalingFactor(vector<Lorentz5Momentum> pout,
			    vector<Energy> mon,Energy roots,
			    Lorentz5Momentum ppartner, Energy mbar,
			    double & k1, double & k2) const {
  ThreeVector<Energy> qtotal;
  vector<Energy2> pmag; 
  for(unsigned int ix=0;ix<pout.size();++ix) {
    pmag.push_back(pout[ix].vect().mag2());
    qtotal+=pout[ix].vect();
  }
  Energy2 dot1 = qtotal*ppartner.vect();
  Energy2 qmag2=qtotal.mag2();
  double a = -dot1/qmag2;
  static const Energy eps=1e-10*GeV;
  unsigned int itry(0);
  Energy numer(ZERO),denom(ZERO);
  k1=1.;
  do {
    ++itry;
    numer=denom=0.*GeV;
    double k12=sqr(k1);
    for(unsigned int ix=0;ix<pout.size();++ix) {
      Energy en = sqrt(pmag[ix]/k12+sqr(mon[ix]));
      numer += en;
      denom += pmag[ix]/en;
    }
    Energy en = sqrt(qmag2/k12+sqr(mbar));
    numer += en-roots;
    denom += qmag2/en;
    k1 += numer/denom*k12*k1;
    if(abs(k1)>1e10) return false;
  }
  while (abs(numer)>eps&&itry<100);
  k1 = abs(k1);
  k2 = a*k1;
  return itry<100;
}

void KinematicsReconstructor::
deconstructInitialFinalSystem(HardTreePtr tree,vector<HardBranchingPtr> jets,
			      ShowerInteraction type) const {
  HardBranchingPtr incoming;
  Lorentz5Momentum pin[2],pout[2],pbeam;
  HardBranchingPtr initial;
  Energy mc(ZERO);
  for(unsigned int ix=0;ix<jets.size();++ix) {
    // final-state parton
    if(jets[ix]->status()==HardBranching::Outgoing) {
      pout[0] += jets[ix]->branchingParticle()->momentum();
      mc = jets[ix]->branchingParticle()->thePEGBase() ? 
	jets[ix]->branchingParticle()->thePEGBase()->mass() :
	jets[ix]->branchingParticle()->dataPtr()->mass();
    }
    // initial-state parton
    else {
      pin[0]  += jets[ix]->branchingParticle()->momentum();
      initial = jets[ix];
      pbeam = jets[ix]->beam()->momentum();
      Energy scale=pbeam.t();
      pbeam = Lorentz5Momentum(ZERO,pbeam.vect().unit()*scale);
      incoming = jets[ix];
      while(incoming->parent()) incoming = incoming->parent();
    }
  }
  if(jets.size()>2) {
    pout[0].rescaleMass();
    mc = pout[0].mass();
  }
  // work out the boost to the Breit frame
  Lorentz5Momentum pa = pout[0]-pin[0];
  Axis axis(pa.vect().unit());
  LorentzRotation rot;
  double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
  if(axis.perp2()>0.) {
    rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    rot.rotateX(Constants::pi);
    rot.boostZ( pa.e()/pa.vect().mag());
  }
  // transverse part
  Lorentz5Momentum paxis=rot*pbeam;
  Boost trans = -1./paxis.e()*paxis.vect();
  trans.setZ(0.);
  rot.boost(trans);
  pa *= rot;
  // reference vectors
  Lorentz5Momentum n1(ZERO,ZERO,-pa.z(),-pa.z());
  Lorentz5Momentum n2(ZERO,ZERO, pa.z(),-pa.z());
  Energy2 n1n2 = n1*n2;
  // decompose the momenta
  Lorentz5Momentum qbp=rot*pin[0],qcp= rot*pout[0];
  double a[2],b[2];
  a[0] = n2*qbp/n1n2;
  b[0] = n1*qbp/n1n2;
  a[1] = n2*qcp/n1n2;
  b[1] = n1*qcp/n1n2;
  Lorentz5Momentum qperp = qbp-a[0]*n1-b[0]*n2;
  // before reshuffling
  Energy Q = abs(pa.z());
  double c = sqr(mc/Q);
  Lorentz5Momentum pb(ZERO,ZERO,0.5*Q*(1.+c),0.5*Q*(1.+c));
  Lorentz5Momentum pc(ZERO,ZERO,0.5*Q*(c-1.),0.5*Q*(1.+c));
  double anew[2],bnew[2];
  anew[0] = pb*n2/n1n2;
  bnew[0] = 0.5*(qbp.m2()-qperp.m2())/n1n2/anew[0];
  bnew[1] = pc*n1/n1n2;
  anew[1] = 0.5*qcp.m2()/bnew[1]/n1n2;
  Lorentz5Momentum qnewb = (anew[0]*n1+bnew[0]*n2+qperp);
  Lorentz5Momentum qnewc = (anew[1]*n1+bnew[1]*n2);
  // initial-state boost
  LorentzRotation rotinv=rot.inverse();
  LorentzRotation transb=rotinv*solveBoostZ(qnewb,qbp)*rot;
  // final-state boost
  LorentzRotation transc=rotinv*solveBoost(qnewc,qcp)*rot;
  // this will need changing for more than one outgoing particle
  // set the pvectors
  for(unsigned int ix=0;ix<jets.size();++ix) {
    if(jets[ix]->status()==HardBranching::Incoming) {
      jets[ix]->pVector(pbeam);
      jets[ix]->showerMomentum(rotinv*pb);
      incoming->pVector(jets[ix]->pVector());
    }
    else {
      jets[ix]->pVector(rotinv*pc);
      jets[ix]->showerMomentum(jets[ix]->pVector());
    }
  }
  // find the colour partners
  ShowerParticleVector particles;
  vector<Lorentz5Momentum> ptemp;
  set<HardBranchingPtr>::const_iterator cjt;
  for(cjt=tree->branchings().begin();cjt!=tree->branchings().end();++cjt) {
    ptemp.push_back((**cjt).branchingParticle()->momentum());
    (**cjt).branchingParticle()->set5Momentum((**cjt).showerMomentum());
    particles.push_back((**cjt).branchingParticle());
  }
  dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->partnerFinder()
    ->setInitialEvolutionScales(particles,false,type,false);
  unsigned int iloc(0);
  for(cjt=tree->branchings().begin();cjt!=tree->branchings().end();++cjt) {
    // reset the momentum
    (**cjt).branchingParticle()->set5Momentum(ptemp[iloc]);
    ++iloc;
  }
  for(vector<HardBranchingPtr>::const_iterator cjt=jets.begin();
      cjt!=jets.end();++cjt) {
    // sort out the partners
    tShowerParticlePtr partner = 
      (*cjt)->branchingParticle()->partner();
    if(!partner) continue;
    tHardBranchingPtr branch;
    for(set<HardBranchingPtr>::const_iterator 
	  clt=tree->branchings().begin();clt!=tree->branchings().end();++clt) {
      if((**clt).branchingParticle()==partner) {
	(**cjt).colourPartner(*clt);
  	branch=*clt;
	break;
      }
    }
    // compute the reference vectors
    // both incoming, should all ready be done
    if((**cjt).status()==HardBranching::Incoming &&
       branch->status()==HardBranching::Incoming) {
      Energy etemp = (*cjt)->beam()->momentum().z();
      Lorentz5Momentum nvect(ZERO, ZERO,-etemp, abs(etemp));
      tHardBranchingPtr branch2 = *cjt;     
      (**cjt).nVector(nvect);
      while (branch2->parent()) {
	branch2=branch2->parent();
	branch2->nVector(nvect);
      }
    }
    // both outgoing
    else if((**cjt).status()==HardBranching::Outgoing&&
	     branch->status()==HardBranching::Outgoing) {
      Boost boost=((*cjt)->pVector()+branch->pVector()).findBoostToCM();
      Lorentz5Momentum pcm = branch->pVector();
      pcm.boost(boost);
      Lorentz5Momentum nvect = Lorentz5Momentum(ZERO,pcm.vect());
      nvect.boost( -boost);
      (**cjt).nVector(nvect);
    }
    else if((**cjt).status()==HardBranching::Incoming) {
      Lorentz5Momentum pa = -(**cjt).showerMomentum()+branch->showerMomentum();
      Lorentz5Momentum pb =  (**cjt).showerMomentum();
      Axis axis(pa.vect().unit());
      LorentzRotation rot;
      double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
      if(axis.perp2()>1e-20) {
	rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
	rot.rotateX(Constants::pi);
      }
      if(abs(1.-pa.e()/pa.vect().mag())>1e-6) rot.boostZ( pa.e()/pa.vect().mag());
      pb*=rot;
      Boost trans = -1./pb.e()*pb.vect();
      trans.setZ(0.);
      rot.boost(trans);
      Energy scale=(**cjt).beam()->momentum().t();
      Lorentz5Momentum pbasis(ZERO,(**cjt).beam()->momentum().vect().unit()*scale);
      Lorentz5Momentum pcm = rot*pbasis;
      rot.invert();
      Lorentz5Momentum nvect = rot*Lorentz5Momentum(ZERO,-pcm.vect());
      (**cjt).nVector(nvect);
      tHardBranchingPtr branch2 = *cjt;     
      while (branch2->parent()) {
	branch2=branch2->parent();
	branch2->nVector(nvect);
      }
    }
    else if(branch->status()==HardBranching::Incoming) {
      Lorentz5Momentum nvect=Lorentz5Momentum(ZERO,branch->showerMomentum().vect());
      (**cjt).nVector(nvect);
    }
  }
  // now compute the new momenta
  for(vector<HardBranchingPtr>::const_iterator cjt=jets.begin();
      cjt!=jets.end();++cjt) {
    if((**cjt).status()==HardBranching::Outgoing) {
      (**cjt).setMomenta(transc,1.,Lorentz5Momentum());
    }
  }
  incoming->setMomenta(transb,1.,Lorentz5Momentum());
}

void KinematicsReconstructor::deepTransform(PPtr particle,
					const LorentzRotation & r,
					bool match,
					PPtr original) const {
  if(_boosts.find(particle)!=_boosts.end()) {
    _boosts[particle].push_back(r);
  }
  Lorentz5Momentum porig = particle->momentum();
  if(!original) original = particle;
  for ( int i = 0, N = particle->children().size(); i < N; ++i ) {
    deepTransform(particle->children()[i],r,
		  particle->children()[i]->id()==original->id()&&match,original);
  }
  particle->transform(r);
  // transform the p and n vectors
  ShowerParticlePtr sparticle = dynamic_ptr_cast<ShowerParticlePtr>(particle);
  if(sparticle && sparticle->showerBasis()) {
    sparticle->showerBasis()->transform(r);
  }
  if ( particle->next() ) deepTransform(particle->next(),r,match,original);
  if(!match) return;
  if(!particle->children().empty()) return;
  // force the mass shell
  if(particle->dataPtr()->stable()) {
    Lorentz5Momentum ptemp = particle->momentum();
    ptemp.rescaleEnergy();
    particle->set5Momentum(ptemp);
  }
  // check if there's a daughter tree which also needs boosting
  map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator tit;
  for(tit  = _currentTree->treelinks().begin();
      tit != _currentTree->treelinks().end();++tit) {
    // if there is, boost it
    if(tit->second.first && tit->second.second==original) {
      Lorentz5Momentum pnew = tit->first->incomingLines().begin()
	->first->progenitor()->momentum();
      pnew *=  tit->first->transform();
      Lorentz5Momentum pdiff = porig-pnew;
      Energy2 test = sqr(pdiff.x()) + sqr(pdiff.y()) + 
   	             sqr(pdiff.z()) + sqr(pdiff.t());
      LorentzRotation rot;
      if(test>1e-6*GeV2) rot = solveBoost(porig,pnew);
      tit->first->transform(r*rot,false);
      _treeBoosts[tit->first].push_back(r*rot);
    }
  }
}

Energy KinematicsReconstructor::findMass(HardBranchingPtr branch) const {
  // KH - 230909 - If the particle has no children then it will 
  // not have showered and so it should be "on-shell" so we can
  // get it's mass from it's momentum. This means that the
  // inverseRescalingFactor doesn't give any nans or do things 
  // it shouldn't if it gets e.g. two Z bosons generated with
  // off-shell masses. This is for sure not the best solution.
  // PR 1/1/10 modification to previous soln
  // PR 28/8/14 change to procedure and factorize into a function
  if(branch->children().empty()) {
    return branch->branchingParticle()->mass();
  }
  else if(!branch->children().empty() &&
	  !branch->branchingParticle()->dataPtr()->stable() ) {
    for(unsigned int ix=0;ix<branch->children().size();++ix) {
      if(branch->branchingParticle()->id()==
	 branch->children()[ix]->branchingParticle()->id())
	return findMass(branch->children()[ix]);
    }
  }
  return branch->branchingParticle()->dataPtr()->mass();
}

vector<double>
KinematicsReconstructor::inverseInitialStateRescaling(double & x1, double & x2,
						  const Lorentz5Momentum & pold,
						  const vector<Lorentz5Momentum> & p,
						  const vector<Lorentz5Momentum> & pq) const {
  // hadronic CMS
  Energy2 s  = (pq[0] +pq[1] ).m2();
  // partonic CMS
  Energy MDY = pold.m();
  // find alpha, beta and pt
  Energy2 p12=pq[0]*pq[1];
  double a[2],b[2];
  Lorentz5Momentum pt[2];
  for(unsigned int ix=0;ix<2;++ix) {
    a[ix] = p[ix]*pq[1]/p12;
    b [ix] = p[ix]*pq[0]/p12;
    pt[ix]    = p[ix]-a[ix]*pq[0]-b[ix]*pq[1];
  }
  // compute kappa
  // we always want to preserve the mass of the system
  double k1(1.),k2(1.);
  if(_initialStateReconOption==0) {
    double rap=pold.rapidity();
    x2 = MDY/sqrt(s*exp(2.*rap));
    x1 = sqr(MDY)/s/x2;
    k1=a[0]/x1;
    k2=b[1]/x2;
  }
  // longitudinal momentum
  else if(_initialStateReconOption==1) {
    double A = 1.;
    double C = -sqr(MDY)/s;
    double B = 2.*pold.z()/sqrt(s);
    if(abs(B)>1e-10) {
      double discrim = 1.-4.*A*C/sqr(B);
      if(discrim < 0.) throw KinematicsReconstructionVeto();
      x1 = B>0. ? 0.5*B/A*(1.+sqrt(discrim)) : 0.5*B/A*(1.-sqrt(discrim));
    }
    else {
      x1 = -C/A;
      if( x1 <= 0.) throw KinematicsReconstructionVeto();
      x1 = sqrt(x1);
    }
    x2 = sqr(MDY)/s/x1;
    k1=a[0]/x1;
    k2=b[1]/x2;
  }
  // preserve mass and don't scale the softer system
  // to reproduce the dipole kinematics
  else if(_initialStateReconOption==2) {
    // in this case kp = k1 or k2 depending on who's the harder guy
    k1 = a[0]*b[1]*s/sqr(MDY);
    if ( pt[0].perp2() < pt[1].perp2() ) swap(k1,k2);
    x1 = a[0]/k1;
    x2 = b[1]/k2;
  }
  else
    assert(false);
  // decompose the momenta
  double anew[2] = {a[0]/k1,a[1]*k2};
  double bnew[2] = {b[0]*k1,b[1]/k2};
  vector<double> boost(2);
  for(unsigned int ix=0;ix<2;++ix) {
    boost[ix] = getBeta(a   [ix]+b   [ix], a[ix]   -b   [ix], 
			anew[ix]+bnew[ix], anew[ix]-bnew[ix]);
  }
  return boost;
}

vector<double>
KinematicsReconstructor::initialStateRescaling(double x1, double x2, 
					   const Lorentz5Momentum & pold,
					   const vector<Lorentz5Momentum> & p,
					   const vector<Lorentz5Momentum> & pq,
					   const vector<Energy>& highestpts) const {
  Energy2 S = (pq[0]+pq[1]).m2();
  // find alphas and betas in terms of desired basis
  Energy2 p12 = pq[0]*pq[1];
  double a[2] = {p[0]*pq[1]/p12,p[1]*pq[1]/p12};
  double b[2] = {p[0]*pq[0]/p12,p[1]*pq[0]/p12};
  Lorentz5Momentum p1p = p[0] - a[0]*pq[0] - b[0]*pq[1];
  Lorentz5Momentum p2p = p[1] - a[1]*pq[0] - b[1]*pq[1];
  // compute kappa
  // we always want to preserve the mass of the system
  Energy MDY = pold.m();
  Energy2 A = a[0]*b[1]*S;
  Energy2 B = Energy2(sqr(MDY)) - (a[0]*b[0]+a[1]*b[1])*S - (p1p+p2p).m2();
  Energy2 C = a[1]*b[0]*S;
  double rad = 1.-4.*A*C/sqr(B);
  if(rad < 0.) throw KinematicsReconstructionVeto();
  double kp = B/(2.*A)*(1.+sqrt(rad));
  // now compute k1
  // conserve rapidity
  double k1(0.);
  double k2(0.);
  if(_initialStateReconOption==0) {
    rad = kp*(b[0]+kp*b[1])/(kp*a[0]+a[1]);
    rad *= pq[0].z()<ZERO ? exp(-2.*pold.rapidity()) : exp(2.*pold.rapidity());
    if(rad <= 0.) throw KinematicsReconstructionVeto();
    k1 = sqrt(rad);
    k2 = kp/k1;
  }
  // conserve longitudinal momentum
  else if(_initialStateReconOption==1) {
    double a2 = (a[0]+a[1]/kp);
    double b2 = -x2+x1;
    double c2 = -(b[1]*kp+b[0]);
    if(abs(b2)>1e-10) {
      double discrim = 1.-4.*a2*c2/sqr(b2);
      if(discrim < 0.) throw KinematicsReconstructionVeto();
      k1 = b2>0. ? 0.5*b2/a2*(1.+sqrt(discrim)) : 0.5*b2/a2*(1.-sqrt(discrim));
    }
    else {
      k1 = -c2/a2;
      if( k1 <= 0.) throw KinematicsReconstructionVeto();
      k1 = sqrt(k1);
    }
    k2 = kp/k1;
  }
  // preserve mass and don't scale the softer system
  // to reproduce the dipole kinematics
  else if(_initialStateReconOption==2) {
    // in this case kp = k1 or k2 depending on who's the harder guy
    k1 = kp; k2 = 1.;
    if ( highestpts[0] < highestpts[1] )
      swap(k1,k2);
  }
  else
    assert(false);
  // calculate the boosts
  vector<double> beta(2);
  beta[0] = getBeta((a[0]+b[0]), (a[0]-b[0]), (k1*a[0]+b[0]/k1), (k1*a[0]-b[0]/k1));
  beta[1] = getBeta((a[1]+b[1]), (a[1]-b[1]), (a[1]/k2+k2*b[1]), (a[1]/k2-k2*b[1]));
  if (pq[0].z() > ZERO) {
    beta[0] = -beta[0];
    beta[1] = -beta[1];
  }
  return beta;
}

void KinematicsReconstructor::
reconstructColourSinglets(vector<ShowerProgenitorPtr> & ShowerHardJets,
			  ShowerInteraction type) const {
  // identify and catagorize the colour singlet systems
  unsigned int nnun(0),nnii(0),nnif(0),nnf(0),nni(0);
  vector<ColourSingletSystem> 
    systems(identifySystems(set<ShowerProgenitorPtr>(ShowerHardJets.begin(),ShowerHardJets.end()),
			    nnun,nnii,nnif,nnf,nni));
  // now decide what to do
  // initial-initial connection and final-state colour singlet systems
  LorentzRotation toRest,fromRest;
  bool applyBoost(false),general(false);
  // Drell-Yan type
  if(nnun==0&&nnii==1&&nnif==0&&nnf>0&&nni==0) {
    // reconstruct initial-initial system
    for(unsigned int ix=0;ix<systems.size();++ix) {
      if(systems[ix].type==II) 
	reconstructInitialInitialSystem(applyBoost,toRest,fromRest,
					systems[ix].jets);
    }
    if(type!=ShowerInteraction::QCD) {
      combineFinalState(systems);
      general=false;
    }
  }
  // DIS and VBF type
  else if(nnun==0&&nnii==0&&((nnif==1&&nnf>0&&nni==1)||
			     (nnif==2&&       nni==0))) {
    // check these systems can be reconstructed
    for(unsigned int ix=0;ix<systems.size();++ix) {
      // compute q^2
      if(systems[ix].type!=IF) continue;
      Lorentz5Momentum q;
      for(unsigned int iy=0;iy<systems[ix].jets.size();++iy) {
	if(systems[ix].jets[iy]->progenitor()->isFinalState())
	  q += systems[ix].jets[iy]->progenitor()->momentum();
	else
	  q -= systems[ix].jets[iy]->progenitor()->momentum();
      }
      q.rescaleMass();
      // check above cut
      if(abs(q.m())>=_minQ) continue;
      if(nnif==1&&nni==1) {
	throw KinematicsReconstructionVeto();
      }
      else {
	general = true;
	break;
      }
    }
    if(!general) {
      for(unsigned int ix=0;ix<systems.size();++ix) {
	if(systems[ix].type==IF)
	  reconstructInitialFinalSystem(systems[ix].jets);
      }
    }
  }
  // e+e- type
  else if(nnun==0&&nnii==0&&nnif==0&&nnf>0&&nni==2) {
    general = type!=ShowerInteraction::QCD;
  }
  // general type
  else {
    general = true;
  }
  // final-state systems except for general recon
  if(!general) {
    for(unsigned int ix=0;ix<systems.size();++ix) {
      if(systems[ix].type==F)
	reconstructFinalStateSystem(applyBoost,toRest,fromRest,
				    systems[ix].jets);
    }
  }
  else {
    reconstructGeneralSystem(ShowerHardJets);
  }
}

void KinematicsReconstructor::findInitialBoost(const Lorentz5Momentum & pold,
					   const Lorentz5Momentum & pnew,
					   LorentzRotation & toRest,
					   LorentzRotation & fromRest) const {
  // do one boost
  if(_initialBoost==0) {
    toRest   = LorentzRotation(pold.findBoostToCM());
    fromRest = LorentzRotation(pnew.boostVector());
  }
  else if(_initialBoost==1) {
    // boost to rest frame
    // first transverse
    toRest = Boost(-pold.x()/pold.t(),-pold.y()/pold.t(),0.);
    // then longitudinal
    double beta = pold.z()/sqrt(pold.m2()+sqr(pold.z()));
    toRest.boost((Boost(0.,0.,-beta)));
    // boost from rest frame
    // first apply longitudinal boost
    beta = pnew.z()/sqrt(pnew.m2()+sqr(pnew.z()));
    fromRest=LorentzRotation(Boost(0.,0.,beta));
    // then transverse one
    fromRest.boost(Boost(pnew.x()/pnew.t(),
			 pnew.y()/pnew.t(),0.));
  }
  else
    assert(false);
}
