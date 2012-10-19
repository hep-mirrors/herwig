// -*- C++ -*-
//
// QTildeReconstructor.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeReconstructor class.
//

#include "QTildeReconstructor.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"
#include <cassert>

using namespace Herwig;

namespace {

enum SystemType { UNDEFINED=-1, II, IF, F ,I };

struct ColourSingletSystem {

  ColourSingletSystem() : type(UNDEFINED) {}

  ColourSingletSystem(SystemType intype,ShowerProgenitorPtr inpart) 
    : type(intype),jets(1,inpart) {}

  /**
   * The type of system
   */
  SystemType type;

  /**
   *  The jets in the system
   */
  vector<ShowerProgenitorPtr> jets;
};

struct ColourSingletShower {

  ColourSingletShower() : type(UNDEFINED) {}

  ColourSingletShower(SystemType intype,HardBranchingPtr inpart) 
    : type(intype),jets(1,inpart) {}

  /**
   * The type of system
   */
  SystemType type;

  /**
   *  The jets in the system
   */
  vector<HardBranchingPtr> jets;

};

}

void QTildeReconstructor::persistentOutput(PersistentOStream & os) const {
  os << _reconopt << _initialBoost << ounit(_minQ,GeV) << _noRescale 
     << _noRescaleVector;
}

void QTildeReconstructor::persistentInput(PersistentIStream & is, int) {
  is >> _reconopt >> _initialBoost >> iunit(_minQ,GeV) >> _noRescale 
     >> _noRescaleVector;  
}

ClassDescription<QTildeReconstructor> QTildeReconstructor::initQTildeReconstructor;
// Definition of the static class description member.

void QTildeReconstructor::Init() {

  static ClassDocumentation<QTildeReconstructor> documentation
    ( "This class is responsible for the kinematics reconstruction of the showering,",
      " including the kinematics reshuffling necessary to compensate for the recoil"
      "of the emissions." );

  static Switch<QTildeReconstructor,unsigned int> interfaceReconstructionOption
    ("ReconstructionOption",
     "Option for the kinematics reconstruction",
     &QTildeReconstructor::_reconopt, 0, false, false);
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

  static Parameter<QTildeReconstructor,Energy> interfaceMinimumQ2
    ("MinimumQ2",
     "The minimum Q2 for the reconstruction of initial-final systems",
     &QTildeReconstructor::_minQ, GeV, 0.001*GeV, 1e-6*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static RefVector<QTildeReconstructor,ParticleData> interfaceNoRescale
    ("NoRescale",
     "Particles which shouldn't be rescaled to be on shell by the shower",
     &QTildeReconstructor::_noRescaleVector, -1, false, false, true, false, false);

  static Switch<QTildeReconstructor,unsigned int> interfaceInitialInitialBoostOption
    ("InitialInitialBoostOption",
     "Option for how the boost from the system before ISR to that after ISR is applied.",
     &QTildeReconstructor::_initialBoost, 0, false, false);
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

}

void QTildeReconstructor::doinit() {
  KinematicsReconstructor::doinit();
  _noRescale = set<cPDPtr>(_noRescaleVector.begin(),_noRescaleVector.end());
}

bool QTildeReconstructor::
reconstructTimeLikeJet(const tShowerParticlePtr particleJetParent,
		       unsigned int iopt) const {
  assert(particleJetParent);
  bool emitted=true;
  // if this is not a fixed point in the reconstruction
  if( !particleJetParent->isReconstructionFixedPoint() ) {
    // if not a reconstruction fixpoint, dig deeper for all children:
    for ( ParticleVector::const_iterator cit = 
	    particleJetParent->children().begin();
	  cit != particleJetParent->children().end(); ++cit )
      reconstructTimeLikeJet(dynamic_ptr_cast<ShowerParticlePtr>(*cit),iopt);
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
	   !_progenitor->data().stable()) {
	  jetGrandParent->showerKinematics()->reconstructLast(particleJetParent,iopt,
							      _progenitor->mass());
	}
	else {
	  jetGrandParent->showerKinematics()->reconstructLast(particleJetParent,iopt);
	}
      }
    }
    // otherwise
    else {
      Energy dm = particleJetParent->data().constituentMass();
      if (abs(dm-particleJetParent->momentum().m())>0.001*MeV
	  &&particleJetParent->dataPtr()->stable()
	  &&particleJetParent->id()!=ParticleID::gamma
	  &&_noRescale.find(particleJetParent->dataPtr())==_noRescale.end()) {
	Lorentz5Momentum dum =  particleJetParent->momentum();
	if(dm>dum.e()) throw KinematicsReconstructionVeto();
	dum.setMass(dm); 
	dum.rescaleRho(); 
	particleJetParent->set5Momentum(dum);  
      } 
      else {
	emitted=false;
      }
    }
  }
  // recursion has reached an endpoint once, ie we can reconstruct the
  // kinematics from the children.
  if( !(particleJetParent->isReconstructionFixedPoint()) ) 
    particleJetParent->showerKinematics()
      ->reconstructParent( particleJetParent, particleJetParent->children() );
  return emitted;
}

bool QTildeReconstructor::
reconstructHardJets(ShowerTreePtr hard,
		    const map<tShowerProgenitorPtr,
		    pair<Energy,double> > & intrinsic) const {
  _currentTree = hard;
  _intrinsic=intrinsic;
  // extract the particles from the ShowerTree
  vector<ShowerProgenitorPtr> ShowerHardJets=hard->extractProgenitors();
  try {
    // old recon method, using new member functions
    if(_reconopt==0) {
      reconstructGeneralSystem(ShowerHardJets);
    }
    // reconstruction based on coloured systems
    else {
      // identify the colour singlet systems
      vector<ColourSingletSystem> systems;
      vector<bool> done(ShowerHardJets.size(),false);
      for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
	// if not treated create new system
	if(done[ix]) continue;
	systems.push_back(ColourSingletSystem(UNDEFINED,ShowerHardJets[ix]));
	done[ix] = true;
	if(!ShowerHardJets[ix]->progenitor()->coloured()) continue;
	// now find the colour connected particles
	vector<unsigned int> iloc(1,ix);
	do {
	  vector<unsigned int> temp=findPartners(iloc.back(),ShowerHardJets);
	  iloc.pop_back();
	  for(unsigned int iy=0;iy<temp.size();++iy) {
	    if(!done[temp[iy]]) {
	      done[temp[iy]] = true;
	      iloc.push_back(temp[iy]);
	      systems.back().jets.push_back(ShowerHardJets[temp[iy]]);
	    }
	  }
	}
	while(!iloc.empty());
      }
      // catagorize the systems
      unsigned int nnun(0),nnii(0),nnif(0),nnf(0),nni(0);
      for(unsigned int ix=0;ix<systems.size();++ix) {
	unsigned int ni(0),nf(0);
	for(unsigned int iy=0;iy<systems[ix].jets.size();++iy) {
	  if(systems[ix].jets[iy]->progenitor()->isFinalState()) ++nf;
	  else                                                   ++ni;
	}
	// type
	// initial-initial
	if(ni==2&&nf==0) {
	  systems[ix].type = II;
	  ++nnii;
	}
	// initial only
	else if(ni==1&&nf==0) {
	  systems[ix].type = I;
	  ++nni;
	}
	// initial-final
	else if(ni==1&&nf>0) {
	  systems[ix].type = IF;
	  ++nnif;
	}
	// final only
	else if(ni==0&&nf>0) {
	  systems[ix].type = F;
	  ++nnf;
	}
	// otherwise unknown
	else {
	  systems[ix].type = UNDEFINED;
	  ++nnun;
	}
      }
      // now decide what to do
      // initial-initial connection and final-state colour singlet systems
      LorentzRotation toRest,fromRest;
      bool applyBoost(false);
      bool general(false);
      // Drell-Yan type
      if(nnun==0&&nnii==1&&nnif==0&&nnf>0&&nni==0) {
	// reconstruct initial-initial system
	for(unsigned int ix=0;ix<systems.size();++ix) {
	  if(systems[ix].type==II) 
	    reconstructInitialInitialSystem(applyBoost,toRest,fromRest,
					    systems[ix].jets);
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
	// only FS needed
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
  }
  catch(KinematicsReconstructionVeto) {
    _progenitor=tShowerParticlePtr();
    _intrinsic.clear();
    _currentTree = tShowerTreePtr();
    return false;
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
      _currentTree = tShowerTreePtr();
      return false;
    }
  }
  _currentTree = tShowerTreePtr();
  return true;
}

double 
QTildeReconstructor::solveKfactor(const Energy & root_s, 
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

bool QTildeReconstructor::
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
      reconstructTimeLikeJet(child,0);
      // calculate the momentum of the particle
      Lorentz5Momentum pnew=p->momentum()-child->momentum();
      pnew.rescaleMass();
      p->children()[0]->set5Momentum(pnew);
    }
  }
  return emitted;
}

Boost QTildeReconstructor::
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

bool QTildeReconstructor::
reconstructDecayJets(ShowerTreePtr decay) const {
  _currentTree = decay;
  try {
    // extract the particles from the ShowerTree
    vector<ShowerProgenitorPtr> ShowerHardJets=decay->extractProgenitors();
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
	nvect = initial->progenitor()->showerKinematics()->getBasis()[1];
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
      atLeastOnce |= reconstructTimeLikeJet(tempJetKin.parent,0);
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
    _currentTree = tShowerTreePtr();
    return false;
  }
  _currentTree = tShowerTreePtr();
  return true;
}

bool QTildeReconstructor::
reconstructDecayJet( const tShowerParticlePtr p) const {
  if(p->children().empty()) return false;
  tShowerParticlePtr child;
  // if branching reconstruct time-like child
  child = dynamic_ptr_cast<ShowerParticlePtr>(p->children()[1]);
  if(child) {
    _progenitor=child;
    reconstructTimeLikeJet(child,1);
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

bool QTildeReconstructor::
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
  Lorentz5Momentum ptest;
  for(unsigned int ix=0;ix<jetKinematics.size();++ix) {
    pmag.push_back(jetKinematics[ix].p.vect().mag2());
    total+=jetKinematics[ix].q.mass();
    ptest+=jetKinematics[ix].p;
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

bool QTildeReconstructor::
deconstructDecayJets(HardTreePtr decay, cEvolverPtr,
		     ShowerInteraction::Type) const {
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
//   // find the evolution partners
//   ShowerParticleVector particles;
//   particles.push_back((**decay->incoming().begin()).branchingParticle());
//   for(cit=branchings.begin();cit!=branchings.end();++cit) {
//     if((**cit).status()==HardBranching::Outgoing)
//       particles.push_back((*cit)->branchingParticle());
//   }
//   // partners should
//   evolver->showerModel()->partnerFinder()
//     ->setInitialEvolutionScales(particles,true,type,false);
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

double QTildeReconstructor::
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
      throw Exception() << "Rescaling factor is imaginary in  QTildeReconstructor::"
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
    vector<Energy> root;
    root.resize(pout.size());
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
  if(isnan(lambda))
    throw Exception() << "Rescaling factor is nan in  QTildeReconstructor::"
		      << "inverseRescalingFactor " 
		      << Exception::eventerror;
  return lambda;
}

bool QTildeReconstructor::
deconstructGeneralSystem(HardTreePtr tree,
			 cEvolverPtr evolver,
			 ShowerInteraction::Type type) const {
  // extract incoming and outgoing particles
  ColourSingletShower in,out;
  for(set<HardBranchingPtr>::const_iterator it=tree->branchings().begin();
      it!=tree->branchings().end();++it) {
    if((**it).status()==HardBranching::Incoming) in .jets.push_back(*it);
    else                  out.jets.push_back(*it);
  }
  // do the initial-state reconstruction
  LorentzRotation toRest,fromRest;
  bool applyBoost(false);
  deconstructInitialInitialSystem(applyBoost,toRest,fromRest,
				  tree,in.jets,type);
  // do the final-state reconstruction
  deconstructFinalStateSystem(toRest,fromRest,tree,
			      out.jets,evolver,type);
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

bool QTildeReconstructor::deconstructHardJets(HardTreePtr tree,
					      cEvolverPtr evolver,
					      ShowerInteraction::Type type) const {
  // inverse of old recon method
  if(_reconopt==0) {
    return deconstructGeneralSystem(tree,evolver,type);
  }
  // inverse of reconstruction based on coloured systems
  else {
    // identify the colour singlet systems
    vector<ColourSingletShower> systems;
    set<HardBranchingPtr> done;
    for(set<HardBranchingPtr>::const_iterator it=tree->branchings().begin();
	it!=tree->branchings().end();++it) {
      // if not treated create new system
      if(done.find(*it)!=done.end()) continue;
      done.insert(*it);
      systems.push_back(ColourSingletShower(UNDEFINED,*it));
      if(!(**it).branchingParticle()->coloured()) continue;
      // now find the colour connected particles
      findPartners(*it,done,tree->branchings(),systems.back().jets);
    }
    // catagorize the systems
    unsigned int nnun(0),nnii(0),nnif(0),nnf(0),nni(0);
    for(unsigned int ix=0;ix<systems.size();++ix) {
      unsigned int ni(0),nf(0);
      for(unsigned int iy=0;iy<systems[ix].jets.size();++iy) {
	if(systems[ix].jets[iy]->status()==HardBranching::Outgoing) ++nf;
	else                                                        ++ni;
      }
      // type
      // initial-initial
      if(ni==2&&nf==0) {
	systems[ix].type = II;
	++nnii;
      }
      // initial only
      else if(ni==1&&nf==0) {
	systems[ix].type = I;
	++nni;
      }
      // initial-final
      else if(ni==1&&nf>0) {
	systems[ix].type = IF;
	++nnif;
      }
      // final only
      else if(ni==0&&nf>0) {
	systems[ix].type = F;
	++nnf;
      }
      // otherwise unknown
      else {
	systems[ix].type = UNDEFINED;
	++nnun;
      }
    }
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
    }
    // DIS and VBF type
    else if(nnun==0&&nnii==0&&((nnif==1&&nnf>0&&nni==1)||
			       (nnif==2&&       nni==0))) {
      for(unsigned int ix=0;ix<systems.size();++ix) {
	if(systems[ix].type==IF)
 	  deconstructInitialFinalSystem(tree,systems[ix].jets,evolver,type);
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
				      systems[ix].jets,evolver,type);
      }
    }
    else {
      return deconstructGeneralSystem(tree,evolver,type);
    }
    return true;
  }
}

vector<unsigned int> QTildeReconstructor::
findPartners(unsigned int iloc ,
	     vector<ShowerProgenitorPtr> jets) const {
  vector<unsigned int> output;
  for(unsigned int iy=0;iy<jets.size();++iy) {
    if(!jets[iy]->progenitor()->data().coloured()||iy==iloc) continue;
    bool isPartner = false;
    // both in either initial or final state
    if(jets[iloc]->progenitor()->isFinalState()!=jets[iy]->progenitor()->isFinalState()) {
      //loop over all the colours of both
      for(unsigned int ix=0; ix<CLSIZE(jets[iloc]); ++ix){
        for(unsigned int jx=0; jx<CLSIZE(jets[iy]); ++jx){
          if(CL(jets[iloc],ix) && CL(jets[iloc],ix)==CL(jets[iy],jx))
            isPartner = true; 
        }
      }    
      if(!isPartner){
        //loop over anti colours of both
        for(unsigned int ix=0; ix<ACLSIZE(jets[iloc]); ++ix){
          for(unsigned int jx=0; jx<ACLSIZE(jets[iy]); ++jx){
            if(ACL(jets[iloc],ix) && ACL(jets[iloc],ix)==ACL(jets[iy],jx))
              isPartner = true; 
          }
        }
      }
    }
    else{
      //loop over the colours of the first and the anti-colours of the other
      for(unsigned int ix=0; ix<CLSIZE(jets[iloc]); ++ix){
        for(unsigned int jx=0; jx<ACLSIZE(jets[iy]); ++jx){
          if(CL(jets[iloc],ix) && CL(jets[iloc],ix)==ACL(jets[iy],jx))
            isPartner = true; 
        }
      }
      if(!isPartner){
        //loop over the anti-colours of the first and the colours of the other
        for(unsigned int ix=0; ix<ACLSIZE(jets[iloc]); ++ix){
          for(unsigned int jx=0; jx<CLSIZE(jets[iy]); jx++){
            if(ACL(jets[iloc],ix) && ACL(jets[iloc],ix)==CL(jets[iy],jx))
              isPartner = true;
          }
        }
      }
    }
    // special for sources/sinks
    if(jets[iloc]->progenitor()->colourLine()) {
      if(jets[iloc]->progenitor()->colourLine()->sourceNeighbours().first) {
        tColinePair lines = jets[iloc]->progenitor()->colourLine()->sourceNeighbours();
        if(lines.first ==  jets[iy]->progenitor()->    colourLine() || 
           lines.first ==  jets[iy]->progenitor()->    colourLine() || 
           lines.second == jets[iy]->progenitor()->antiColourLine() || 
           lines.second == jets[iy]->progenitor()->antiColourLine()) 
          isPartner = true;
      }
      if(jets[iloc]->progenitor()->colourLine()->sinkNeighbours().first) {
        tColinePair lines = jets[iloc]->progenitor()->colourLine()->sinkNeighbours();
        if(lines.first  == jets[iy]->progenitor()->    colourLine() || 
           lines.first  == jets[iy]->progenitor()->    colourLine() || 
           lines.second == jets[iy]->progenitor()->antiColourLine() || 
           lines.second == jets[iy]->progenitor()->antiColourLine())
          isPartner = true;
      }
    }
    if(jets[iloc]->progenitor()->antiColourLine()) {
      if(jets[iloc]->progenitor()->antiColourLine()->sourceNeighbours().first) {
        tColinePair lines = jets[iloc]->progenitor()->antiColourLine()->sourceNeighbours();
        if(lines.first  == jets[iy]->progenitor()->    colourLine() || 
           lines.first  == jets[iy]->progenitor()->    colourLine() || 
           lines.second == jets[iy]->progenitor()->antiColourLine() || 
           lines.second == jets[iy]->progenitor()->antiColourLine())
          isPartner = true;
      }
      if(jets[iloc]->progenitor()->antiColourLine()->sinkNeighbours().first) {
        tColinePair lines = jets[iloc]->progenitor()->antiColourLine()->sinkNeighbours();
        if(lines.first  == jets[iy]->progenitor()->    colourLine() || 
           lines.first  == jets[iy]->progenitor()->    colourLine() || 
           lines.second == jets[iy]->progenitor()->antiColourLine() || 
           lines.second == jets[iy]->progenitor()->antiColourLine())
          isPartner = true;
      }
    }
    if(isPartner)
      output.push_back(iy);
  }
  return output;
}

void QTildeReconstructor::
reconstructInitialFinalSystem(vector<ShowerProgenitorPtr> jets) const {
  Lorentz5Momentum pin[2],pout[2];
  bool atLeastOnce(false);
  for(unsigned int ix=0;ix<jets.size();++ix) {
    // final-state parton
    if(jets[ix]->progenitor()->isFinalState()) {
      pout[0] +=jets[ix]->progenitor()->momentum();
      _progenitor = jets[ix]->progenitor();
      atLeastOnce |= reconstructTimeLikeJet(jets[ix]->progenitor(),0);
    }
    // initial-state parton
    else {
      pin[0]  +=jets[ix]->progenitor()->momentum();
      atLeastOnce |= reconstructSpaceLikeJet(jets[ix]->progenitor());
      assert(!jets[ix]->original()->parents().empty());
    }
  }
  // add intrinsic pt if needed
  atLeastOnce |= addIntrinsicPt(jets);
  // momenta after showering
  for(unsigned int ix=0;ix<jets.size();++ix) {
    if(jets[ix]->progenitor()->isFinalState())
      pout[1] += jets[ix]->progenitor()->momentum();
    else
      pin[1]  += jets[ix]->progenitor()->momentum();
  }
  // work out the boost to the Breit frame
  Lorentz5Momentum pa = pout[0]-pin[0];
  Lorentz5Momentum pb = pin[0];
  Axis axis(pa.vect().unit());
  LorentzRotation rot;
  double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
  rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  rot.rotateX(Constants::pi);
  rot.boostZ( pa.e()/pa.vect().mag());
  Lorentz5Momentum ptemp=rot*pb;
  Boost trans = -1./ptemp.e()*ptemp.vect();
  trans.setZ(0.);
  rot.boost(trans);
  pa *=rot;
  // project and calculate rescaling
  // reference vectors
  Lorentz5Momentum n1(ZERO,ZERO,-pa.z(),-pa.z());
  Lorentz5Momentum n2(ZERO,ZERO, pa.z(),-pa.z());
  Energy2 n1n2 = n1*n2;
  // decompose the momenta
  Lorentz5Momentum qbp=rot*pin[1],qcp= rot*pout[1];
  double a[2],b[2];
  a[0] = n2*qbp/n1n2;
  b[0] = n1*qbp/n1n2;
  Lorentz5Momentum qperp = qbp-a[0]*n1-b[0]*n2;
  a[1] = 0.5*(qcp.m2()-qperp.m2())/n1n2;
  b[1] = 1.;
  double A(0.5*a[0]),B(b[0]*a[0]-a[1]*b[1]-0.25),C(-0.5*b[0]);
  if(sqr(B)-4.*A*C<0.) throw KinematicsReconstructionVeto();
  double kb = 0.5*(-B+sqrt(sqr(B)-4.*A*C))/A;
  double kc = (a[0]*kb-0.5)/a[1];
  if(kc==0.) throw KinematicsReconstructionVeto();
  Lorentz5Momentum pnew[2] = { a[0]*kb*n1+b[0]/kb*n2+qperp,
			       a[1]*kc*n1+b[1]/kc*n2+qperp};
  LorentzRotation rotinv=rot.inverse();
  LorentzRotation transb=rotinv*solveBoostZ(pnew[0],qbp)*rot;
  LorentzRotation transc=rotinv*solveBoost(pnew[1],qcp)*rot;
  for(unsigned int ix=0;ix<jets.size();++ix) {
    if(jets[ix]->progenitor()->isFinalState())
      deepTransform(jets[ix]->progenitor(),transc);
    else {
      tPPtr parent;
      boostChain(jets[ix]->progenitor(),transb,parent);
    }
  }
}

bool QTildeReconstructor::addIntrinsicPt(vector<ShowerProgenitorPtr> jets) const {
  bool added=false;
  // add the intrinsic pt if needed
  for(unsigned int ix=0;ix<jets.size();++ix) {
    // only for initial-state particles which haven't radiated
    if(jets[ix]->progenitor()->isFinalState()||
       jets[ix]->hasEmitted()) continue;
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

LorentzRotation QTildeReconstructor::
solveBoost(const double k, const Lorentz5Momentum & newq, 
	   const Lorentz5Momentum & oldp ) const {
  Energy q = newq.vect().mag(); 
  Energy2 qs = sqr(q); 
  Energy2 Q2 = newq.mass2(); 
  Energy kp = k*(oldp.vect().mag()); 
  Energy2 kps = sqr(kp);
  double betam = (q*newq.e() - kp*sqrt(kps + Q2))/(kps + qs + Q2); 
  Boost beta = -betam*(k/kp)*oldp.vect();
  // note that (k/kp)*oldp.vect() = oldp.vect()/oldp.vect().mag() but cheaper. 
  ThreeVector<Energy2> ax = newq.vect().cross( oldp.vect() ); 
  double delta = newq.vect().angle( oldp.vect() );
  LorentzRotation R;
  using Constants::pi;
  if ( ax.mag2()/GeV2/MeV2 > 1e-16 ) {
    R.rotate( delta, unitVector(ax) ).boost( beta );
  } 
  else if(abs(delta-pi)/pi < 0.001) {
    double phi=2.*pi*UseRandom::rnd();
    Axis axis(cos(phi),sin(phi),0.);
    axis.rotateUz(newq.vect().unit());
    R.rotate(delta,axis).boost( beta );
  }
  else {
    R.boost( beta );
  } 
  return R;
}

LorentzRotation QTildeReconstructor::solveBoost(const Lorentz5Momentum & q, 
						const Lorentz5Momentum & p ) const {
  Energy modp = p.vect().mag();
  Energy modq = q.vect().mag();
  double betam = (p.e()*modp-q.e()*modq)/(sqr(modq)+sqr(modp)+p.mass2());
  Boost beta = -betam*q.vect().unit();
  ThreeVector<Energy2> ax = p.vect().cross( q.vect() ); 
  double delta = p.vect().angle( q.vect() );
  LorentzRotation R;
  using Constants::pi;
  if ( ax.mag2()/GeV2/MeV2 > 1e-16 ) {
    R.rotate( delta, unitVector(ax) ).boost( beta );
  } 
  else {
    R.boost( beta );
  } 
  return R;
}

LorentzRotation QTildeReconstructor::solveBoostZ(const Lorentz5Momentum & q, 
						 const Lorentz5Momentum & p ) const {
  static const Energy2 eps2 = 1e-8*GeV2;
  static const Energy  eps  = 1e-4 *GeV;
  LorentzRotation R;
  double beta;
  Energy2 den = (p.t()*q.t()-p.z()*q.z());
  Energy2 num = -(p.z()*q.t()-q.z()*p.t());
  if(abs(den)<eps2||abs(num)<eps2) {
    if(abs(p.t()-abs(p.z()))<eps&&abs(q.t()-abs(q.z()))<eps) {
      double ratio = sqr(q.t()/p.t());
      beta  = -(1.-ratio)/(1.+ratio); 
    }
    else {
      beta=0.;
    }
  }
  else {
    beta = num/den;
    
  }
  R.boostZ(beta);
  Lorentz5Momentum ptest = R*p;
  if(ptest.z()/q.z() < 0. || ptest.t()/q.t() < 0. ) {
    throw KinematicsReconstructionVeto();
  }
  return R;
}

void QTildeReconstructor::
reconstructFinalStateSystem(bool applyBoost, 
			    const LorentzRotation &   toRest,
			    const LorentzRotation & fromRest, 
			    vector<ShowerProgenitorPtr> jets) const {
  // special for case of individual particle
  if(jets.size()==1) {
    LorentzRotation trans(toRest);
    trans.transform(fromRest);
    deepTransform(jets[0]->progenitor(),trans);
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
  // check if in CMF frame
  Boost beta_cm = pcm.findBoostToCM();
  bool gottaBoost = (beta_cm.mag() > 1e-12);
  // collection of pointers to initial hard particle and jet momenta
  // for final boosts
  JetKinVect jetKinematics;
  vector<ShowerProgenitorPtr>::const_iterator cit;
  for(cit = jets.begin(); cit != jets.end(); cit++) {
    JetKinStruct tempJetKin;      
    tempJetKin.parent = (*cit)->progenitor(); 
    if(gottaBoost) {
      tempJetKin.parent->boost(beta_cm);
      map<tShowerTreePtr,pair<tShowerProgenitorPtr,
	tShowerParticlePtr> >::const_iterator tit;
      for(tit  = _currentTree->treelinks().begin();
	  tit != _currentTree->treelinks().end();++tit) {
	if(tit->second.first && tit->second.second==tempJetKin.parent)
	  tit->first->transform(LorentzRotation(beta_cm),false);
      }
    }
    tempJetKin.p = (*cit)->progenitor()->momentum();
    _progenitor=tempJetKin.parent;
    radiated |= reconstructTimeLikeJet((*cit)->progenitor(),0);
    tempJetKin.q = (*cit)->progenitor()->momentum();
    jetKinematics.push_back(tempJetKin);
  }
  // find the rescaling factor
  double k = 0.0;
  if(radiated) {
    k = solveKfactor(pcm.m(), jetKinematics);
  }
  // perform the rescaling and boosts
  for(JetKinVect::iterator it = jetKinematics.begin();
      it != jetKinematics.end(); ++it) {
    LorentzRotation Trafo = LorentzRotation(); 
    if(radiated) Trafo = solveBoost(k, it->q, it->p);
    if(gottaBoost) Trafo.boost(-beta_cm);
    if(applyBoost) {
      Trafo.transform(  toRest);
      Trafo.transform(fromRest);
    }
    if(radiated || gottaBoost || applyBoost) deepTransform(it->parent,Trafo);
  }
}

void QTildeReconstructor::
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
  }
  // check if intrinsic pt to be added
  radiated |= !_intrinsic.empty();
  // if no radiation return
  if(!radiated) return;
  // initial state shuffling
  applyBoost=false;
  vector<Lorentz5Momentum> p, pq, p_in;
  for(unsigned int ix=0;ix<jets.size();++ix) {
    // at momentum to vector
    p_in.push_back(jets[ix]->progenitor()->momentum());
    // reconstruct the jet
    radiated |= reconstructSpaceLikeJet(jets[ix]->progenitor());
    assert(!jets[ix]->original()->parents().empty());
    Energy etemp = jets[ix]->original()->parents()[0]->momentum().z();
    Lorentz5Momentum ptemp = Lorentz5Momentum(ZERO, ZERO, etemp, abs(etemp));
    pq.push_back(ptemp);
  }
  // add the intrinsic pt if needed
  radiated |=addIntrinsicPt(jets);
  for(unsigned int ix=0;ix<jets.size();++ix) {
    p.push_back(jets[ix]->progenitor()->momentum());
  }
  double x1 = p_in[0].z()/pq[0].z();
  double x2 = p_in[1].z()/pq[1].z();
  Energy MDY = (p_in[0] + p_in[1]).m();
  Energy2 S = (pq[0]+pq[1]).m2();
  // if not need don't apply boosts
  if(!(radiated && p.size() == 2 && pq.size() == 2)) return;
  applyBoost=true;
  // find alphas and betas in terms of desired basis      
  Energy2 p12 = pq[0]*pq[1];
  double a[2] = {p[0]*pq[1]/p12,p[1]*pq[1]/p12};
  double b[2] = {p[0]*pq[0]/p12,p[1]*pq[0]/p12};
  Lorentz5Momentum p1p = p[0] - a[0]*pq[0] - b[0]*pq[1];
  Lorentz5Momentum p2p = p[1] - a[1]*pq[0] - b[1]*pq[1];
  // compute kappa
  Energy2 A = a[0]*b[1]*S;
  Energy2 B = Energy2(sqr(MDY)) - (a[0]*b[0]+a[1]*b[1])*S - (p1p+p2p).m2();
  Energy2 C = a[1]*b[0]*S; 
  double rad = 1.-4.*A*C/sqr(B);
  if(rad < 0.) throw KinematicsReconstructionVeto();
  double kp = B/(2.*A)*(1.+sqrt(rad));
  // now compute k1, k2
  rad = kp*(b[0]+kp*b[1])/(kp*a[0]+a[1])*(x1/x2);  
  if(rad <= 0.) throw KinematicsReconstructionVeto();
  double k1 = sqrt(rad);
  double k2 = kp/k1;
  double beta[2] = 
    {getBeta((a[0]+b[0]), (a[0]-b[0]), (k1*a[0]+b[0]/k1), (k1*a[0]-b[0]/k1)),
     getBeta((a[1]+b[1]), (a[1]-b[1]), (a[1]/k2+k2*b[1]), (a[1]/k2-k2*b[1]))};
  if (pq[0].z() > ZERO) {
    beta[0] = -beta[0]; 
    beta[1] = -beta[1];
  }
  // apply the boosts
  Lorentz5Momentum newcmf;
  for(unsigned int ix=0;ix<jets.size();++ix) {
    tPPtr toBoost = jets[ix]->progenitor();
    Boost betaboost(0, 0, beta[ix]);
    tPPtr parent;
    boostChain(toBoost, LorentzRotation(betaboost),parent);
    if(parent->momentum().e()/pq[ix].e()>1.||
       parent->momentum().z()/pq[ix].z()>1.) throw KinematicsReconstructionVeto();
    newcmf+=toBoost->momentum();
  }
  if(newcmf.m()<ZERO||newcmf.e()<ZERO) throw KinematicsReconstructionVeto();
  // do one boost
  toRest   = LorentzRotation(pcm.findBoostToCM());
  if(_initialBoost==0) {
    fromRest = LorentzRotation(newcmf.boostVector());
  }
  else if(_initialBoost==1) {
    // first apply longitudinal boost
    double beta = newcmf.z()/sqrt(newcmf.m2()+sqr(newcmf.z()));
    fromRest=LorentzRotation(Boost(0.,0.,beta));
    // then transverse one
    Energy pT = sqrt(sqr(newcmf.x())+sqr(newcmf.y()));
    beta = pT/newcmf.t();
    fromRest.boost(Boost(beta*newcmf.x()/pT,beta*newcmf.y()/pT,0.));
  }
  else
    assert(false);
}

void QTildeReconstructor::
deconstructInitialInitialSystem(bool & applyBoost,
				LorentzRotation & toRest,
				LorentzRotation & fromRest,
				HardTreePtr tree,
				vector<HardBranchingPtr> jets,
				ShowerInteraction::Type) const {
  // get the momenta of the particles
  vector<Lorentz5Momentum> pin;
  vector<Lorentz5Momentum> pq;
  vector<HardBranchingPtr>::iterator cit;
  for(cit=jets.begin();cit!=jets.end();++cit) {
    pin.push_back((*cit)->branchingParticle()->momentum());
    Energy etemp = (*cit)->beam()->momentum().z();
    pq.push_back(Lorentz5Momentum(ZERO, ZERO,etemp, abs(etemp)));
  }
  bool order = (*tree->incoming().begin())->beam()->momentum().z()/pq[0].z()<0.;
  assert(pin.size()==2);
  // decompose the momenta
  double alpha[2],beta[2];
  Energy2 p12=pq[0]*pq[1];
  Lorentz5Momentum pt[2];
  for(unsigned int ix=0;ix<2;++ix) {
    alpha[ix] = pin[ix]*pq[1]/p12;
    beta [ix] = pin[ix]*pq[0]/p12;
    pt[ix]    = pin[ix]-alpha[ix]*pq[0]-beta[ix]*pq[1];
  }
  // parton level centre-of-mass
  Lorentz5Momentum pcm=pin[0]+pin[1];
  pcm.rescaleMass();
  double rap=pcm.rapidity();
  // hadron level cmf
  Energy2 s  = (pq[0] +pq[1] ).m2();
  // calculate the x values 
  double x0 = sqrt(pcm.mass2()/s*exp(2.*rap));
  double x[2]={x0, pcm.mass2()/s/x0};
  if(pq[0].z()<ZERO) swap(x[0],x[1]);
  double k1=alpha[0]/x[0],k2=beta[1]/x[1];
  double alphanew[2]={alpha[0]/k1,alpha[1]*k2};
  double betanew [2]={beta [0]*k1,beta [1]/k2};
  double boost[2];
  for(unsigned int ix=0;ix<2;++ix) {
    boost[ix] = getBeta(alpha   [ix]+beta   [ix], alpha[ix]   -beta   [ix], 
			alphanew[ix]+betanew[ix], alphanew[ix]-betanew[ix]);
    if (pq[0].z() > ZERO) beta[ix]*=-1.;
  }
  // apply the boost the the particles
  // first incoming particle
  if(order) swap(pq[0],pq[1]);
  // now apply the boosts
  Boost betaboost(0.,0.,-boost[0]);
  LorentzRotation R;
  R.boost(betaboost);
  set<HardBranchingPtr>::const_iterator cjt=tree->incoming().begin();
  (*cjt)->pVector(pq[0]);
  (*cjt)->nVector(pq[1]);
  (*cjt)->setMomenta(R,1.,Lorentz5Momentum());
  // second incoming particle
  betaboost = Boost(0.,0.,-boost[1]);
  R=LorentzRotation(betaboost);
  ++cjt;
  (*cjt)->pVector(pq[1]);
  (*cjt)->nVector(pq[0]);
  (*cjt)->setMomenta(R,1.,Lorentz5Momentum());
  jets[0]->showerMomentum(x[0]*jets[0]->pVector());
  jets[1]->showerMomentum(x[1]*jets[1]->pVector());
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

void QTildeReconstructor::
deconstructFinalStateSystem(const LorentzRotation &   toRest,
			    const LorentzRotation & fromRest,
			    HardTreePtr tree, vector<HardBranchingPtr> jets,
			    cEvolverPtr evolver,
			    ShowerInteraction::Type type) const {
  if(jets.size()==1) {
    LorentzRotation R(toRest);
    R.transform(fromRest);
    // \todo What does this do?    tree->showerRot( R );
    jets[0]->original(R*jets[0]->branchingParticle()->momentum());
    jets[0]->showerMomentum(R*jets[0]->branchingParticle()->momentum());
    // find the colour partners
    ShowerParticleVector particles;
    vector<Lorentz5Momentum> ptemp;
    set<HardBranchingPtr>::const_iterator cjt;
    for(cjt=tree->branchings().begin();cjt!=tree->branchings().end();++cjt) {
      ptemp.push_back((**cjt).branchingParticle()->momentum());
      (**cjt).branchingParticle()->set5Momentum((**cjt).showerMomentum());
      particles.push_back((**cjt).branchingParticle());
    }
    evolver->showerModel()->partnerFinder()
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
  for(cit=jets.begin();cit!=jets.end();++cit) {
    pout.push_back((*cit)->branchingParticle()->momentum());
    // KH - 230909 - If the particle has no children then it will 
    // not have showered and so it should be "on-shell" so we can
    // get it's mass from it's momentum. This means that the
    // inverseRescalingFactor doesn't give any nans or do things 
    // it shouldn't if it gets e.g. two Z bosons generated with
    // off-shell masses. This is for sure not the best solution.
    // PR 1/1/10 modification to previous soln
    if((*cit)->branchingParticle()->children().size()==0 ||
       (!(*cit)->branchingParticle()->dataPtr()->coloured() &&
	!(*cit)->branchingParticle()->dataPtr()->stable()) ) 
	mon.push_back(pout.back().mass());
    else
      mon.push_back((*cit)->branchingParticle()->dataPtr()->mass());
  }
  // boost all the momenta to the rest frame of the decaying particle
  Lorentz5Momentum pin;
  for(unsigned int ix=0;ix<pout.size();++ix) {
    pout[ix].transform(toRest);
    pin += pout[ix];
  }
  pin.rescaleMass();
  // rescaling factor
  double lambda=inverseRescalingFactor(pout,mon,pin.mass());
  // now calculate the p reference vectors 
  for(cit=jets.begin();cit!=jets.end();++cit){
    Lorentz5Momentum pvect = (*cit)->branchingParticle()->momentum();
    pvect.transform(toRest);
    pvect /= lambda;
    if((*cit)->branchingParticle()->children().size()==0 ||
       (!(*cit)->branchingParticle()->dataPtr()->coloured() &&
	!(*cit)->branchingParticle()->dataPtr()->stable()) ) 
      pvect.setMass((*cit)->branchingParticle()->momentum().mass());
    else
      pvect.setMass((*cit)->branchingParticle()->dataPtr()->mass());
    pvect.rescaleEnergy();
    pvect.transform(fromRest);
    (*cit)->pVector(pvect);
    (*cit)->showerMomentum(pvect);
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
  evolver->showerModel()->partnerFinder()
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
    else if(!(**cjt).status()==HardBranching::Incoming&&
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

Energy QTildeReconstructor::momConsEq(const double & k, 
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

void QTildeReconstructor::boostChain(tPPtr p, const LorentzRotation &bv,
				     tPPtr & parent) const {
  if(!p->parents().empty()) boostChain(p->parents()[0], bv,parent);
  else parent=p;
  p->transform(bv);
  if(p->children().size()==2) {
    if(dynamic_ptr_cast<ShowerParticlePtr>(p->children()[1]))
      deepTransform(p->children()[1],bv);
  }
}

void QTildeReconstructor::
reconstructGeneralSystem(vector<ShowerProgenitorPtr> & ShowerHardJets) const {
  // general recon, all initial-state in one system and final-state
  // in another
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
  reconstructInitialInitialSystem(applyBoost,toRest,fromRest,in.jets);
  // reconstruct the final-state systems
  reconstructFinalStateSystem(applyBoost,toRest,fromRest,out.jets);
}

bool QTildeReconstructor::
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

void QTildeReconstructor::
findPartners(HardBranchingPtr branch,set<HardBranchingPtr> & done,
             const set<HardBranchingPtr> & branchings,
             vector<HardBranchingPtr> & jets) const {
  tShowerParticlePtr part=branch->branchingParticle();
  unsigned int partNumColourLines  = 
    branch->branchingParticle()->colourInfo()->    colourLines().size();
  unsigned int partNumAColourLines =
    branch->branchingParticle()->colourInfo()->antiColourLines().size();
  for(set<HardBranchingPtr>::const_iterator cit=branchings.begin();
      cit!=branchings.end();++cit) {
    if(done.find(*cit)!=done.end()||!(**cit).branchingParticle()->coloured())
      continue;
    bool isPartner = false;
    // one initial and one final
    if(branch->status()!=(**cit).status()) {
      if(part->colourLine()) {
        for(unsigned int ix=0; ix<partNumColourLines; ++ix){
          for(unsigned int jx=0; jx<CLSIZE(cit); ++jx){
            if(part->colourInfo()->colourLines()[ix] == 
               (**cit).branchingParticle()->colourInfo()->colourLines()[jx]){
              isPartner = true;
              break;
            }
          }
        }
      }
      if(part->antiColourLine()) {
        for(unsigned int ix=0; ix<partNumAColourLines; ++ix){
          for(unsigned int jx=0; jx<ACLSIZE(cit); ++jx){
            if(part->colourInfo()->antiColourLines()[ix] == 
               (**cit).branchingParticle()->colourInfo()->antiColourLines()[jx]){
              isPartner = true;
              break;
            }
          }
        }
      }
    }
    // both in either initial or final state
    else {
      if(part->colourLine()) {
        for(unsigned int ix=0; ix<partNumColourLines; ++ix) {
          for(unsigned int jx=0; jx<ACLSIZE(cit); ++jx) {
            if(part->colourInfo()->colourLines()[ix] == 
               (**cit).branchingParticle()->colourInfo()->antiColourLines()[jx]){
              isPartner = true;
              break;
            }
          }
        }
      }
      if(part->antiColourLine()) {
        for(unsigned int ix=0; ix<partNumAColourLines; ++ix){
          for(unsigned int jx=0; jx<CLSIZE(cit); ++jx){
            if(part->colourInfo()->antiColourLines()[ix] == 
               (**cit).branchingParticle()->colourInfo()->colourLines()[jx]){
            isPartner = true;
            break;
            }
          }
        }
      }
    }
    if(isPartner) {
      jets.push_back(*cit);
      done.insert(*cit);
      findPartners(*cit,done,branchings,jets);
    }
  }
}

void QTildeReconstructor::
deconstructInitialFinalSystem(HardTreePtr tree,vector<HardBranchingPtr> jets,
			      cEvolverPtr evolver,
			      ShowerInteraction::Type type) const {
  HardBranchingPtr incoming;
  Lorentz5Momentum pin[2],pout[2],pbeam;
  HardBranchingPtr initial;
  Energy mc(ZERO);
  for(unsigned int ix=0;ix<jets.size();++ix) {
    // final-state parton
    if(jets[ix]->status()==HardBranching::Outgoing) {
      pout[0] += jets[ix]->branchingParticle()->momentum();
      mc = jets[ix]->branchingParticle()->getThePEGBase() ? 
	jets[ix]->branchingParticle()->getThePEGBase()->mass() :
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
  evolver->showerModel()->partnerFinder()
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


void QTildeReconstructor::deepTransform(PPtr particle,
					const LorentzRotation & r,
					bool match,
					PPtr original) const {
  Lorentz5Momentum porig = particle->momentum();
  if(!original) original = particle;
  for ( int i = 0, N = particle->children().size(); i < N; ++i ) {
    deepTransform(particle->children()[i],r,
		  particle->children()[i]->id()==original->id()&&match,original);
  }
  particle->transform(r);
  if ( particle->next() ) deepTransform(particle->next(),r,match,original);
  if(!match) return;
  if(!particle->children().empty()) return;
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
    }
  }
}


