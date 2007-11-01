// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeReconstructor class.
//

#include "QTildeReconstructor.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/Timer.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"

using namespace Herwig;

NoPIOClassDescription<QTildeReconstructor> QTildeReconstructor::initQTildeReconstructor;
// Definition of the static class description member.

void QTildeReconstructor::Init() {

  static ClassDocumentation<QTildeReconstructor> documentation
    ( "This class is responsible for the kinematics reconstruction of the showering,",
      " including the kinematics reshuffling necessary to compensate for the recoil"
      "of the emissions." );

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
  Timer<1100> timer("QTildeReconstructor::reconstructHardJets");
  try {
    bool radiated[2] = {false,false};
    // find the hard process centre-of-mass energy
    Lorentz5Momentum p_cm[2] = {Lorentz5Momentum(),Lorentz5Momentum()};
    // create a vector of the hard particles
    map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator mit;
    vector<ShowerProgenitorPtr> ShowerHardJets;
    for(mit=hard->incomingLines().begin();mit!=hard->incomingLines().end();++mit)
      ShowerHardJets.push_back((*mit).first);
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mjt;
    for(mjt=hard->outgoingLines().begin();mjt!=hard->outgoingLines().end();++mjt)
      ShowerHardJets.push_back((*mjt).first);
    for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
      // final-state jet
      if(ShowerHardJets[ix]->progenitor()->isFinalState()) {
	// did it radiate
	radiated[1] |=ShowerHardJets[ix]->hasEmitted();
	// add momentum
	p_cm[1]+=ShowerHardJets[ix]->progenitor()->momentum();
      }
      // initial-state jet
      else {
	// did it radiate
	radiated[0]|=ShowerHardJets[ix]->hasEmitted();
	// add momentum
	p_cm[0]+=ShowerHardJets[ix]->progenitor()->getThePEGBase()->momentum();
      }
    }
    radiated[0]|=!intrinsic.empty();
    // initial state shuffling
    // the boosts for the initial state
    Boost boostRest,boostNewF;
    bool applyBoost(false);
    if(radiated[0])
      applyBoost=reconstructISJets(p_cm[0],ShowerHardJets,intrinsic,
				   boostRest,boostNewF);
    // final-state reconstruction
    // check if in CMF frame
    Boost beta_cm = p_cm[1].findBoostToCM();
    bool gottaBoost = (beta_cm.mag() > 1e-12);
    // check if any radiation
    bool atLeastOnce = radiated[1];
    // collection of pointers to initial hard particle and jet momenta
    // for final boosts
    JetKinVect jetKinematics;
    vector<ShowerProgenitorPtr>::const_iterator cit;
    for(cit = ShowerHardJets.begin(); cit != ShowerHardJets.end(); cit++) {
      if(!(*cit)->progenitor()->isFinalState()) continue;
      JetKinStruct tempJetKin;      
      tempJetKin.parent = (*cit)->progenitor(); 
      if(gottaBoost) tempJetKin.parent->boost(beta_cm); 
      tempJetKin.p = (*cit)->progenitor()->momentum();
      _progenitor=tempJetKin.parent;
      atLeastOnce |= reconstructTimeLikeJet((*cit)->progenitor(),0);
      tempJetKin.q = (*cit)->progenitor()->momentum();
      jetKinematics.push_back(tempJetKin);
    }
    // find the rescaling factor
    double k = 0.0;
    if(atLeastOnce) {
      k = solveKfactor(p_cm[1].mag(), jetKinematics);
      if(k< 0.) return false;
    }
    // perform the rescaling and boosts
    for(JetKinVect::iterator it = jetKinematics.begin();
	it != jetKinematics.end(); ++it) {
      LorentzRotation Trafo = LorentzRotation(); 
      if(atLeastOnce) Trafo = solveBoost(k, it->q, it->p);
      if(gottaBoost) Trafo.boost(-beta_cm);
      if(atLeastOnce || gottaBoost) it->parent->deepTransform(Trafo);
      if(applyBoost) {
	it->parent->deepBoost(boostRest);
	it->parent->deepBoost(boostNewF);
      }
    }
  }
  catch(KinematicsReconstructionVeto) {
    return false;
  }
  return true;
}

const double 
QTildeReconstructor::solveKfactor(const Energy & root_s, 
				      const JetKinVect & jets) const
{
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
reconstructISJets(Lorentz5Momentum pcm,
		  const vector<ShowerProgenitorPtr> & ShowerHardJets,
		  map<tShowerProgenitorPtr,pair<Energy,double> > intrinsic,
		  Boost & boostRest, Boost & boostNewF) const {
  bool atLeastOnce = false;
  vector<Lorentz5Momentum> p, pq, p_in;
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    // only look at initial state particles
    if(ShowerHardJets[ix]->progenitor()->isFinalState()) continue;
    // at momentum to vector
    p_in.push_back(ShowerHardJets[ix]->progenitor()->getThePEGBase()->momentum());
    // reconstruct the jet
    atLeastOnce |= reconstructSpaceLikeJet(ShowerHardJets[ix]->progenitor());
    p.push_back(ShowerHardJets[ix]->progenitor()->momentum());
    assert(!ShowerHardJets[ix]->original()->parents().empty());
    Energy etemp = ShowerHardJets[ix]->original()->
      parents()[0]->momentum().z();
    Lorentz5Momentum ptemp = Lorentz5Momentum(0*MeV, 0*MeV, etemp, abs(etemp));
    pq.push_back(ptemp);
  }
  // add the intrinsic pt if needed
  int iloc=-1;
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    //dont do anything for the moment for secondary scatters
    if( !ShowerHandler::currentHandler()->FirstInt() ) break;
    // only for initial-state particles which haven't radiated
    if(ShowerHardJets[ix]->progenitor()->isFinalState()) continue;
    ++iloc;
    if(ShowerHardJets[ix]->hasEmitted()) continue;
    if(intrinsic.find(ShowerHardJets[ix])==intrinsic.end()) continue;
    pair<Energy,double> pt=intrinsic[ShowerHardJets[ix]];
    Lorentz5Momentum p_basis(pq[0]),n_basis(pq[1]);
    if(iloc==1) swap(p_basis,n_basis);
    double alpha = ShowerHardJets[ix]->progenitor()->x();
    double beta  = 0.5*(sqr(ShowerHardJets[ix]->progenitor()->data().mass())+
			sqr(pt.first))/alpha/(p_basis*n_basis);
    Lorentz5Momentum pnew=alpha*p_basis+beta*n_basis;
    pnew.setX(pt.first*cos(pt.second));
    pnew.setY(pt.first*sin(pt.second));
    pnew.rescaleMass();
    p[iloc]=pnew;
    ShowerHardJets[ix]->progenitor()->set5Momentum(pnew);
    atLeastOnce=true;
  }
  double x1 = p_in[0].z()/pq[0].z();
  double x2 = p_in[1].z()/pq[1].z();
  Energy MDY = (p_in[0] + p_in[1]).m();
  Energy2 S = (pq[0]+pq[1]).m2();
  // if not need don't apply boosts
  if(!(atLeastOnce && p.size() == 2 && pq.size() == 2)) return false;
  // find alphas and betas in terms of desired basis      
  Energy2 p12 = pq[0]*pq[1];
  double a1 = p[0]*pq[1]/p12;
  double b1 = p[0]*pq[0]/p12;
  double a2 = p[1]*pq[1]/p12;
  double b2 = p[1]*pq[0]/p12;
  Lorentz5Momentum p1p = p[0] - a1*pq[0] - b1*pq[1];
  Lorentz5Momentum p2p = p[1] - a2*pq[0] - b2*pq[1];
  // compute kappa12
  // DGRELL is this textbook method for solving a quadratic
  // numerically stable if 4AC ~= B^2 ? check Numerical Recipes
  double kp = 1.0;
  Energy2 A = a1*b2*S;
  Energy2 B = Energy2(sqr(MDY)) - (a1*b1+a2*b2)*S - (p1p+p2p).mag2();
  Energy2 C = a2*b1*S; 
  double rad = 1.-4.*A*C/sqr(B);
  if (rad >= 0) {kp = B/(2.*A)*(1.+sqrt(rad));}
  else {
    throw Exception() << "QTildeReconstructor::reconstructISJets " 
		      << "WARNING! Can't get kappa_pm!\n"
		      << Exception::eventerror;
  }
  // now compute k1, k2
  double k1 = 1.0, k2 = 1.0;
  rad = kp*(b1+kp*b2)/(kp*a1+a2)*(x1/x2);   
  if (rad > 0) {
    k1 = sqrt(rad);
    k2 = kp/k1;
  } 
  else
    throw Exception() << "QTildeReconstructor::reconstructISJets " 
		      << "  Plus:  k1 = " << k1 
		      << "WARNING! Can't get k1p, k2p!\n"
		      << Exception::eventerror;
  double beta1 = getBeta((a1+b1), (a1-b1), 
			 (k1*a1+b1/k1), (k1*a1-b1/k1));
  double beta2 = getBeta((a2+b2), (a2-b2), 
			 (a2/k2+k2*b2), (a2/k2-k2*b2));
  if (pq[0].z() > 0*MeV) {
    beta1 = -beta1; 
    beta2 = -beta2;
  }
  tPVector toBoost;
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    if(!ShowerHardJets[ix]->progenitor()->isFinalState())
      toBoost.push_back(ShowerHardJets[ix]->progenitor());
  }
  // before boost
  Boost betaboost(0, 0, beta1);
  tPPtr parent;
  boostChain(toBoost[0], betaboost,parent);
  if(parent->momentum().e()/pq[0].e()>1.||parent->momentum().z()/pq[0].z()>1.) throw KinematicsReconstructionVeto();
  betaboost = Boost(0, 0, beta2);
  boostChain(toBoost[1], betaboost,parent);
  if(parent->momentum().e()/pq[1].e()>1.||parent->momentum().z()/pq[1].z()>1.) throw KinematicsReconstructionVeto();
  boostRest = pcm.findBoostToCM();
  Lorentz5Momentum newcmf=(toBoost[0]->momentum() + toBoost[1]->momentum());
  if(newcmf.m()<0.*GeV) throw KinematicsReconstructionVeto();
  boostNewF = newcmf.boostVector();
  return true;
}

bool QTildeReconstructor::
reconstructDecayJets(ShowerTreePtr decay) const {
  try {
    // extract the particles from the ShowerTree
    map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator mit;
    vector<ShowerProgenitorPtr> ShowerHardJets;
    for(mit=decay->incomingLines().begin();mit!=decay->incomingLines().end();++mit)
      ShowerHardJets.push_back((*mit).first);
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mjt;
    for(mjt=decay->outgoingLines().begin();mjt!=decay->outgoingLines().end();++mjt)
      ShowerHardJets.push_back((*mjt).first);
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
    for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
      // only consider final-state jets
      if(!ShowerHardJets[ix]->progenitor()->isFinalState()) continue;
      // do the reconstruction
      JetKinStruct tempJetKin;      
      tempJetKin.parent = ShowerHardJets[ix]->progenitor();
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
