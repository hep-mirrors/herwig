// -*- C++ -*-
//
// QTildeReconstructor.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
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
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"

using namespace Herwig;

namespace {

  enum SystemType { UNDEFINED=-1, II, IF, F ,I };

struct ColourSingletSystem {

  ColourSingletSystem() : type(UNDEFINED) {};

  ColourSingletSystem(SystemType intype,ShowerProgenitorPtr inpart) 
    : type(intype),jets(1,inpart) {};

  /**
   * The type of system
   */
  SystemType type;

  /**
   *  The jets in the system
   */
  vector<ShowerProgenitorPtr> jets;
};

}

void QTildeReconstructor::persistentOutput(PersistentOStream & os) const {
  os << _reconopt;
}

void QTildeReconstructor::persistentInput(PersistentIStream & is, int) {
  is >> _reconopt;  
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

}

bool QTildeReconstructor::
reconstructTimeLikeJet(const tShowerParticlePtr particleJetParent,
		       unsigned int iopt) const {
  if(!particleJetParent)
    throw Exception() << "must have a particle in Kinematics"
		      << "Reconstructor::reconstructTimeLikeJet"
		      << Exception::eventerror;
  bool emitted=true;
  // if this is not a fixed point in the reconstruction
  if( !particleJetParent->isReconstructionFixedPoint() ) {
    // if not a reconstruction fixpoint, dig deeper for all children:
    for ( ParticleVector::const_iterator cit = particleJetParent->children().begin();
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
	  //if (abs(dm-particleJetParent->momentum().mass())>0.05*MeV
	  &&particleJetParent->dataPtr()->stable()
	  &&particleJetParent->id()!=ParticleID::gamma) {
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
		    map<tShowerProgenitorPtr,pair<Energy,double> > intrinsic) const {
  _intrinsic=intrinsic;
  // extract the particles from the ShowerTree
  vector<ShowerProgenitorPtr> ShowerHardJets=hard->extractProgenitors();
  try {
    // old recon method, using new member functions
    if(_reconopt==0) {
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
      Boost toRest,fromRest;
      bool applyBoost(false);
      reconstructInitialInitialSystem(applyBoost,toRest,fromRest,in.jets);
      // reconstruct the final-state systems
      reconstructFinalStateSystem(applyBoost,toRest,fromRest,out.jets);
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
	if(!ShowerHardJets[ix]->progenitor()->coloured()) continue;
	// now find the colour connected particles
	done[ix] = true;
	vector<unsigned int> iloc(1,ix);
	do {
	  vector<unsigned int> temp=findPartners(iloc.back(),ShowerHardJets);
	  done[iloc.back()] = true;
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
      // Drell-Yan type
      if(nnun==0&&nnii==1&&nnif==0&&nnf>0&&nni==0) {
	// reconstruct initial-initial system
	Boost toRest,fromRest;
	bool applyBoost(false);
	for(unsigned int ix=0;ix<systems.size();++ix) {
	  if(systems[ix].type==II) 
	    reconstructInitialInitialSystem(applyBoost,toRest,fromRest,systems[ix].jets);
	}
	// reconstruct the final-state systems
	if(nnf>1) {
	  for(unsigned int ix=0;ix<systems.size();++ix) {
	    if(systems[ix].type==F) 
	      reconstructFinalStateSystem(applyBoost,toRest,fromRest,systems[ix].jets);
	  }
	}
	else {
	  for(unsigned int ix=0;ix<systems.size();++ix) {
	    if(systems[ix].type==F) 
	      reconstructFinalStateSystem(applyBoost,toRest,fromRest,systems[ix].jets);
	  }
	}
      }
      // DIS type
      else if(nnun==0&&nnii==0&&nnif==1&&nnf>0&&nni==1) {
	for(unsigned int ix=0;ix<systems.size();++ix) {
	  if(systems[ix].type==IF) 
	    reconstructInitialFinalSystem(systems[ix].jets);
	}
      }
      // e+e- type
      else if(nnun==0&&nnii==0&&nnif==0&&nnf>0&&nni==2) {
	Boost toRest,fromRest;
	bool applyBoost(false);
	for(unsigned int ix=0;ix<systems.size();++ix) {
	  if(systems[ix].type==F) 
	    reconstructFinalStateSystem(applyBoost,toRest,fromRest,systems[ix].jets);
	}
      }
      // VBF type
      else if(nnun==0&&nnii==0&&nnif==2&&nni==0) {
	// initial-final systems
	for(unsigned int ix=0;ix<systems.size();++ix) {
	  if(systems[ix].type==IF) 
	    reconstructInitialFinalSystem(systems[ix].jets);
	}
	// final-state systems
	Boost toRest,fromRest;
	bool applyBoost(false);
	for(unsigned int ix=0;ix<systems.size();++ix) {
	  if(systems[ix].type==F) 
	    reconstructFinalStateSystem(applyBoost,toRest,fromRest,systems[ix].jets);
	}
      }
      // general type
      else {
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
	Boost toRest,fromRest;
	bool applyBoost(false);
	reconstructInitialInitialSystem(applyBoost,toRest,fromRest,in.jets);
	// reconstruct the final-state systems
	reconstructFinalStateSystem(applyBoost,toRest,fromRest,out.jets);
      }
    }
  }
  catch(KinematicsReconstructionVeto) {
    return false;
  }
  _progenitor=tShowerParticlePtr();
  _intrinsic.clear();
  return true;
}

double 
QTildeReconstructor::solveKfactor(const Energy & root_s, 
				  const JetKinVect & jets) const {
  Energy2 s = sqr(root_s);
  // must be at least two jets
  if ( jets.size() < 2) return -1.0;
  // sum of jet masses must be less than roots
  if(momConsEq( 0.0, root_s, jets )>0.0*MeV) return -1.0;
  // if two jets simple solution
  if ( jets.size() == 2 ) {
    if ( sqr((jets[0].p.x()+jets[1].p.x())/MeV) < 1.e-4 &&
	 sqr((jets[0].p.y()+jets[1].p.y())/MeV) < 1.e-4 &&
	 sqr((jets[0].p.z()+jets[1].p.z())/MeV) < 1.e-4 ) {
      return sqrt( ( sqr(s - jets[0].q.m2() - jets[1].q.m2()) 
		     - 4.*jets[0].q.m2()*jets[1].q.m2() )
		   /(4.*s*jets[0].p.vect().mag2()) );
    } 
    else return -1;
  }
  // i.e. jets.size() > 2, numerically
  // check convergence, if it's a problem maybe use Newton iteration?
  else {
    double k1 = 0.,k2 = 1.,k = 0.; 
    
    if ( momConsEq( k1, root_s, jets ) < 0.0*MeV ) {
      while ( momConsEq( k2, root_s, jets ) < 0.0*MeV ) {
	k1 = k2; 
	k2 *= 2;       
      }
      while ( fabs( (k1 - k2)/(k1 + k2) ) > 1.e-10 ) { 
	if( momConsEq( k2, root_s, jets ) == 0.*MeV ) {
	  return k2; 
	} else {
	  k = (k1+k2)/2.;
	  if ( momConsEq( k, root_s, jets ) > 0*MeV ) {
	    k2 = k;
	  } else {
	    k1 = k; 
	  } 
	}
      }
      return k1; 	  
    } else return -1.;
  }
  return -1.; 
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
solveBoostBeta( const double k, const Lorentz5Momentum & newq, const Lorentz5Momentum & oldp ) {

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
    // if initial state radiation reconsturct the jet and set up the basis vectors
    bool initialrad=false;
    Lorentz5Momentum pjet;
    Lorentz5Momentum nvect;
    ShowerParticlePtr partner;
    Lorentz5Momentum ppartner[2];
    if(radiated[0]) {
      // find the partner
      partner=initial->progenitor()->
	partners()[initial->progenitor()->showerKinematics()->
		   splittingFn()->interactionType()];
      if(partner) ppartner[0]=partner->momentum();
      // reconstruct the decay jet
      initialrad=true;
      reconstructDecayJet(initial->progenitor());
      // momentum of decaying particle after ISR
      pjet=initial->progenitor()->momentum()
	-decay->incomingLines().begin()->second->momentum();
      pjet.rescaleMass();
      // get the n reference vector
      nvect= initial->progenitor()->showerKinematics()->getBasis()[1];
    }
    // find boost to the rest frame if needed
    Boost boosttorest=-initial->progenitor()->momentum().boostVector();
    double gammarest =
      initial->progenitor()->momentum().e()/
      initial->progenitor()->momentum().mass();
    // check if need to boost to rest frame
    bool gottaBoost = (boosttorest.mag() > 1e-12);
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
    Energy inmass(0.*GeV);
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
      if(gottaBoost) tempJetKin.parent->deepTransform(restboost); 
      tempJetKin.q = ShowerHardJets[ix]->progenitor()->momentum();
      jetKinematics.push_back(tempJetKin);
      // check if potential partner of the decay particle
      ShowerParticlePtr ptemp=ShowerHardJets[ix]->progenitor()->partners()
	[ShowerIndex::QCD];
      if(ptemp&&!partner&&!ptemp->isFinalState()) 
	possiblepartners.push_back(tempJetKin);
    }
    // now select the partner of the decaying particle if needed
    if(!partner&&!possiblepartners.empty()) {
      unsigned int iloc = UseRandom::irnd(0,possiblepartners.size()-1);
      partner = possiblepartners[iloc].parent;
      nvect = possiblepartners[iloc].p;
      nvect = Lorentz5Momentum(0.*MeV,0.5*initial->progenitor()->mass()*
			       nvect.vect().unit());
      nvect.boost(-boosttorest,gammarest);
      ppartner[0] = possiblepartners[iloc].p;
    }
    if(partner) ppartner[1]=partner->momentum();
    // calculate the rescaling parameters
    double k1,k2;
    Lorentz5Momentum qt;
    if(!solveDecayKFactor(initial->progenitor()->mass(),nvect,pjet,
			  jetKinematics,partner,ppartner,k1,k2,qt)) return false;
    // apply boosts and rescalings to final-state jets
    for(JetKinVect::iterator it = jetKinematics.begin(); 
	it != jetKinematics.end(); ++it) {
      LorentzRotation Trafo = LorentzRotation(); 
      if(it->parent!=partner) {
	// boost for rescaling
	if(atLeastOnce) {
	  if(it->parent->children().empty()) {
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
	if(atLeastOnce || gottaBoost) it->parent->deepTransform(Trafo);
      }
      else {
	Lorentz5Momentum pnew=ppartner[0];
	pnew *=k1;
	pnew-=qt;
	pnew.setMass(ppartner[1].mass());
	pnew.rescaleEnergy();
	LorentzRotation Trafo=solveBoost(1.,ppartner[1],pnew);
	if(gottaBoost) Trafo.boost(-boosttorest,gammarest);
	partner->deepTransform(Trafo);
      }
    }
  }
  catch(KinematicsReconstructionVeto) {
    return false;
  } 
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
solveDecayKFactor(Energy mb, Lorentz5Momentum n, Lorentz5Momentum pjet, 
		  const JetKinVect & jetKinematics, ShowerParticlePtr partner, 
		  Lorentz5Momentum ppartner[2],
		  double & k1, double & k2,Lorentz5Momentum & qt) const {
  Energy2 pjn  = partner ? pjet.vect()*n.vect()        : 0.*MeV2;
  Energy2 pcn  = partner ? ppartner[0].vect()*n.vect() : 1.*MeV2;
  Energy2 nmag = n.vect().mag2();
  Lorentz5Momentum pn=(pjn/nmag)*n;
  qt=pjet-pn;qt.setE(0.*MeV);
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
    ds    = 0.*MeV;
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
    ++ix;
  }
  while(abs(mb-roots)>eps&&ix<100);
  k1=d1;
  k2=d2;
  // return true if N-R succeed, otherwise false
  return ix<100;
}

bool QTildeReconstructor::reconstructDecayShower(NasonTreePtr decay,
						 EvolverPtr evolver) const {
  // extract the momenta of the particles
  vector<Lorentz5Momentum> pin;
  vector<Lorentz5Momentum> pout;
  vector<Energy> mon;
  set<NasonBranchingPtr>::iterator cit;
  set<NasonBranchingPtr> branchings=decay->branchings();
  // The for-loop goes over all _progenitors_ and stores their momenta
  // and _on-shell_ masses (for quarks it's the current mass), so e.g.
  // for 2 body decays it has 3 loops (1 for each decay product and 1 
  // for the decaying particle). The momenta need not be on-shell. For 
  // input from Nason hardest emission generators one progenitor should 
  // be off-shell and the others should be on-shell.
  for(cit=branchings.begin();cit!=branchings.end();++cit){
    if((*cit)->_particle->isFinalState()) {
      pout.push_back((*cit)->_particle->momentum());
      mon.push_back((*cit)->_particle->dataPtr()->mass());
    }
    else {
      pin.push_back((*cit)->_particle->momentum());
    }
  }
  assert(pin.size()==1);
  // boost all the momenta to the rest frame of the decaying particle
  Boost boostv=-pin[0].boostVector();
  for(unsigned int ix=0;ix<pout.size();++ix) pout[ix].boost(boostv);
  double lambda;
  if(pout.size()==2) { 
    Energy rts(pin[0].m()); 
    double mu_q1(pout[0].m()/rts), mu_q2(pout[1].m()/rts);
    double mu_p1(mon[0]/rts)     , mu_p2(mon[1]/rts);
    lambda = sqrt(
       ((1.+mu_q1+mu_q2)*(1.-mu_q1-mu_q2)*(mu_q1-1.-mu_q2)*(mu_q2-1.-mu_q1))
      /((1.+mu_p1+mu_p2)*(1.-mu_p1-mu_p2)*(mu_p1-1.-mu_p2)*(mu_p2-1.-mu_p1)));
  } else {
    lambda=inverseRescaleingFactor(pout,mon,pin[0].mass());
  }
  if(isnan(lambda)) {
    cerr << "\n\n\nQTildeReconstructor::reconstructDecayShower \n";
    cerr << "lambda = " << lambda << "\n";
    cerr << "particles in the branchings including any children:\n";
    for(cit=branchings.begin();cit!=branchings.end();++cit) {
      cerr << "(*cit)->_particle\n" << *((*cit)->_particle) << "\n";
      if((*cit)->_children.size()!=0) {
	cerr << "Has " << (*cit)->_children.size() << " children:\n";
	for(unsigned  int ix=0;ix< (*cit)->_children.size();++ix) {
	  cerr << *((*cit)->_children[ix]->_particle) << "\n";
	}
      } else { cerr << "No children.\n"; }
    }
    throw Exception() << "recon fails " << Exception::eventerror;
  }
  // now calculate the p reference vectors 
  for(cit=branchings.begin();cit!=branchings.end();++cit){
    if(!(*cit)->_particle->isFinalState()) continue;
    (*cit)->_p=(*cit)->_particle->momentum();
    (*cit)->_p.boost(boostv);
    (*cit)->_p/=lambda;
    (*cit)->_p.setMass((*cit)->_particle->dataPtr()->mass());
    (*cit)->_p.rescaleEnergy();
    (*cit)->_shower = (*cit)->_p;
    (*cit)->_shower.boost(-boostv);
  }
  // find the colour partners
  ShowerParticleVector particles;
  for(cit=branchings.begin();cit!=branchings.end();++cit) {
    particles.push_back((*cit)->_particle);
  }
  evolver->showerModel()->partnerFinder()->setQCDInitialEvolutionScales(particles,true);
  // calculate the reference vectors
  for(cit=branchings.begin();cit!=branchings.end();++cit){
    // find the partner branchings
    tShowerParticlePtr partner=(*cit)->_particle->partners()[ShowerIndex::QCD];
    if(!partner) continue;
    tNasonBranchingPtr branch;
    set<NasonBranchingPtr>::iterator cjt;
    for(cjt=branchings.begin();cjt!=branchings.end();++cjt){
      if(cjt==cit) continue;
      if((*cjt)->_particle==partner) {
 	branch=*cjt;
 	break;
      }
    }
    // If there are only two final-state particles this boost should do
    // nothing, since then we should already have (*cit)->_p=-branch->_p.
    Boost boost=((*cit)->_p+branch->_p).findBoostToCM();
    Lorentz5Momentum pcm = branch->_p;
    pcm.boost(boost);
    (*cit)->_n = Lorentz5Momentum(0.*MeV,pcm.vect());
    (*cit)->_n.boost( -boost);
  }
  // now compute the new momenta 
  for(cit=branchings.begin();cit!=branchings.end();++cit){
    if(!(*cit)->_particle->isFinalState()) continue;
    Energy2 dot=(*cit)->_p*(*cit)->_n;
    double beta = 0.5*((*cit)->_particle->momentum().m2()-sqr((*cit)->_p.mass()))/dot;
    Lorentz5Momentum qnew=(*cit)->_p+beta*(*cit)->_n;
    qnew.rescaleMass();
    // qnew is the unshuffled momentum in the rest frame of the _p basis vectors,
    // for the simple case Z->q qbar g this was checked against analytic formulae.
    // compute the boost
    LorentzRotation A=LorentzRotation(boostv);
    LorentzRotation R=solveBoost(qnew,A*(*cit)->_particle->momentum())*A;
    // when R is applied to (*cit)->_particle->momentum() you get qnew (checked).
    (*cit)->setMomenta(R,1.0,Lorentz5Momentum());
  }
  return true;
}

double QTildeReconstructor::
inverseRescaleingFactor(vector<Lorentz5Momentum> pout,
			vector<Energy> mon,Energy roots) const {
  unsigned int ntry=0;
  double lambda=1.;
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
    Energy sum(0.*MeV);
    for(unsigned int ix=0;ix<pout.size();++ix) {
      root[ix] = sqrt(pmag[ix]/sqr(lambda)+sqr(mon[ix]));
      sum+=root[ix];
    }
    // if accuracy reached exit
    if(abs(sum/roots-1.)<1e-10) break;
    // use Newton-Raphson to compute new guess for lambda
    Energy numer(0.*MeV),denom(0.*MeV);
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
  return lambda;
}


bool QTildeReconstructor::reconstructHardShower(NasonTreePtr hard,EvolverPtr) const {
  // extract the incoming particles
  vector<Lorentz5Momentum> pin;
  vector<Lorentz5Momentum> pout;
  vector<Lorentz5Momentum> pq;
  set<NasonBranchingPtr>::iterator cit;
  set<NasonBranchingPtr> branchings=hard->branchings();
  for(cit=branchings.begin();cit!=branchings.end();++cit){
    if((*cit)->_particle->isFinalState()) {
      pout.push_back((*cit)->_particle->momentum());
    }
    else {
      pin.push_back((*cit)->_particle->momentum());
      Energy etemp = (*cit)->_beam->momentum().z();
      pq.push_back(Lorentz5Momentum(0*MeV, 
				    0*MeV, 
				    etemp, 
				    abs(etemp)));
    }
  }
  bool order = (*hard->incoming().begin())->_beam->momentum().z()/pq[0].z()<0.;
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
  double x[2]={sqrt(pcm.mass2()/s*exp(2.*rap)),pcm.mass2()/s/x[0]};
  if(pq[0].z()<0*MeV) swap(x[0],x[1]);
  double k1=alpha[0]/x[0],k2=beta[1]/x[1];
  double alphanew[2]={alpha[0]/k1,alpha[1]*k2};
  double betanew [2]={beta [0]*k1,beta [1]/k2};
  double boost[2];
  for(unsigned int ix=0;ix<2;++ix) {
    boost[ix] = getBeta(alpha   [ix]+beta   [ix], alpha[ix]   -beta   [ix], 
			alphanew[ix]+betanew[ix], alphanew[ix]-betanew[ix]);
    if (pq[0].z() > 0*MeV) beta[ix]*=-1.;
  }
  // apply the boost the the particles
  // first incoming particle
  if(order) {
    swap(pq[0],pq[1]);
  }
  // now apply the boosts
  Boost betaboost(0.,0.,boost[0]);
  LorentzRotation R;
  R.boost(betaboost);
  branchings=hard->incoming();
  cit=branchings.begin();
  (*cit)->_p=pq[0];
  (*cit)->_n=pq[1];
  (*cit)->setMomenta(R,1.,Lorentz5Momentum());
  // second incoming particle
  betaboost = Boost(0.,0.,boost[1]);
  R=LorentzRotation(betaboost);
  ++cit;
  (*cit)->_p=pq[1];
  (*cit)->_n=pq[0];
  (*cit)->setMomenta(R,1.,Lorentz5Momentum());
  return true;
}

vector<unsigned int>  QTildeReconstructor::findPartners(unsigned int iloc ,
							vector<ShowerProgenitorPtr> jets) const {
  vector<unsigned int> output;
  for(unsigned int iy=0;iy<jets.size();++iy) {
    if(!jets[iy]->progenitor()->data().coloured()||iy==iloc) continue;
    bool isPartner = false;
    // both in either initial or final state
    if(jets[iloc]->progenitor()->isFinalState()!=jets[iy]->progenitor()->isFinalState()) {
      if(jets[iloc]->progenitor()->colourLine() &&
	 jets[iloc]->progenitor()->colourLine() == jets[iy]->progenitor()->colourLine())
	isPartner = true;
      if(jets[iloc]->progenitor()->antiColourLine() &&
	 jets[iloc]->progenitor()->antiColourLine() == jets[iy]->progenitor()->antiColourLine())
	isPartner = true;
    }
    else {
      if(jets[iloc]->progenitor()->colourLine() &&
	 jets[iloc]->progenitor()->colourLine() == jets[iy]->progenitor()->antiColourLine())
	isPartner = true;
      if(jets[iloc]->progenitor()->antiColourLine() &&
	 jets[iloc]->progenitor()->antiColourLine() == jets[iy]->progenitor()->colourLine())
	isPartner = true;
    }
    if(isPartner) output.push_back(iy);
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
  Lorentz5Momentum pc = pout[0];
  Axis axis(pa.vect().unit());
  LorentzRotation rot;
  double sinth(sqrt(1.-sqr(axis.z())));
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
  Lorentz5Momentum n1(0.*MeV,0.*MeV,-pa.z(),-pa.z());
  Lorentz5Momentum n2(0.*MeV,0.*MeV, pa.z(),-pa.z());
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
  Lorentz5Momentum pnew[2] = { a[0]*kb*n1+b[0]/kb*n2+qperp,
			       a[1]*kc*n1+b[1]/kc*n2+qperp};
  LorentzRotation rotinv=rot.inverse();
  LorentzRotation transb=rotinv*solveBoost(pnew[0],qbp)*rot;
  LorentzRotation transc=rotinv*solveBoost(pnew[1],qcp)*rot;
  for(unsigned int ix=0;ix<jets.size();++ix) {
    if(jets[ix]->progenitor()->isFinalState())
      jets[ix]->progenitor()->deepTransform(transc);
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
      p_basis(0*MeV, 0*MeV, etemp, abs(etemp)),
      n_basis(0*MeV, 0*MeV,-etemp, abs(etemp));
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

LorentzRotation QTildeReconstructor::solveBoost(const Lorentz5Momentum & q, 
						const Lorentz5Momentum & p ) const {
  Energy modp = p.vect().mag();
  Energy modq = q.vect().mag();
  double betam = (p.e()*modp-q.e()*modq)/(sqr(modq)+sqr(modp)+p.mass2());
  Boost beta = -betam*q.vect().unit();
  Vector3<Energy2> ax = p.vect().cross( q.vect() ); 
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

void QTildeReconstructor::
reconstructFinalStateSystem(bool applyBoost, Boost toRest, Boost fromRest, 
			    vector<ShowerProgenitorPtr> jets) const {
  // special for case of individual particle
  if(jets.size()==1) {
    jets[0]->progenitor()->deepBoost(toRest);
    jets[0]->progenitor()->deepBoost(fromRest);
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
    if(gottaBoost) tempJetKin.parent->boost(beta_cm); 
    tempJetKin.p = (*cit)->progenitor()->momentum();
    _progenitor=tempJetKin.parent;
    radiated |= reconstructTimeLikeJet((*cit)->progenitor(),0);
    tempJetKin.q = (*cit)->progenitor()->momentum();
    jetKinematics.push_back(tempJetKin);
  }
  // find the rescaling factor
  double k = 0.0;
  if(radiated) {
    k = solveKfactor(pcm.mag(), jetKinematics);
    if(k< 0.) throw KinematicsReconstructionVeto();
  }
  // perform the rescaling and boosts
  for(JetKinVect::iterator it = jetKinematics.begin();
      it != jetKinematics.end(); ++it) {
    LorentzRotation Trafo = LorentzRotation(); 
    if(radiated) Trafo = solveBoost(k, it->q, it->p);
    if(gottaBoost) Trafo.boost(-beta_cm);
    if(radiated || gottaBoost) it->parent->deepTransform(Trafo);
    if(applyBoost) {
      it->parent->deepBoost(toRest);
      it->parent->deepBoost(fromRest);
    }
  }
}

void QTildeReconstructor::
reconstructInitialInitialSystem(bool & applyBoost, Boost & toRest, Boost & fromRest,  
				vector<ShowerProgenitorPtr> jets) const {
  bool radiated = false;
  Lorentz5Momentum pcm;
  // check whether particles radiated and calculate total momentum
  for(unsigned int ix=0;ix<jets.size();++ix) {
    radiated |= jets[ix]->hasEmitted();
    pcm += jets[ix]->progenitor()->getThePEGBase()->momentum();
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
    p_in.push_back(jets[ix]->progenitor()->getThePEGBase()->momentum());
    // reconstruct the jet
    radiated |= reconstructSpaceLikeJet(jets[ix]->progenitor());
    assert(!jets[ix]->original()->parents().empty());
    Energy etemp = jets[ix]->original()->parents()[0]->momentum().z();
    Lorentz5Momentum ptemp = Lorentz5Momentum(0*MeV, 0*MeV, etemp, abs(etemp));
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
  // compute kappa[0]2
  // DGRELL is this textbook method for solving a quadratic
  // numerically stable if 4AC ~= B^2 ? check Numerical Recipes
  Energy2 A = a[0]*b[1]*S;
  Energy2 B = Energy2(sqr(MDY)) - (a[0]*b[0]+a[1]*b[1])*S - (p1p+p2p).mag2();
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
  if (pq[0].z() > 0*MeV) {
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
  if(newcmf.m()<0.*GeV||newcmf.e()<0.*GeV) throw KinematicsReconstructionVeto();
  toRest   = pcm.findBoostToCM();
  fromRest = newcmf.boostVector();
}
