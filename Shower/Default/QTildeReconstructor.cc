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
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
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
    {throw Exception() << "must have a particle in Kinematics"
		       << "Reconstructor::reconstructTimeLikeJet"
		       << Exception::eventerror;}
  bool emitted=true;
  // if this is not a fixed point in the reconstruction
  if( !particleJetParent->isReconstructionFixedPoint() ) {
    // if not a reconstruction fixpoint, dig deeper for all children:
    for ( ParticleVector::const_iterator cit = particleJetParent->children().begin();
	  cit != particleJetParent->children().end(); ++cit )
      {reconstructTimeLikeJet(dynamic_ptr_cast<ShowerParticlePtr>(*cit),iopt);}
  }
  // it is a reconstruction fixpoint, ie kinematical data has to be available 
  else {
    // check if the parent was part of the shower
    ShowerParticlePtr jetGrandParent;
    if(!particleJetParent->parents().empty())
      jetGrandParent= dynamic_ptr_cast<ShowerParticlePtr>
	(particleJetParent->parents()[0]);
    // update if so
    if (jetGrandParent) 
      {if (jetGrandParent->showerKinematics())
	  jetGrandParent->showerKinematics()->reconstructLast(particleJetParent,iopt);}
    // otherwise
    else {
      Energy dm = particleJetParent->data().constituentMass();
      if (abs(dm-particleJetParent->momentum().m())>0.001*MeV
	  //if (abs(dm-particleJetParent->momentum().mass())>0.05*MeV
	  &&particleJetParent->dataPtr()->stable()
	  &&particleJetParent->id()!=ParticleID::gamma) {
	Lorentz5Momentum dum =  particleJetParent->momentum();
	if(dm>dum.e()) throw Veto();
	dum.setMass(dm); 
	dum.rescaleRho(); 
	particleJetParent->set5Momentum(dum);  
      } 
      else {emitted=false;}
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
    Vector3 boostRest,boostNewF;
    bool applyBoost(false);
    if(radiated[0])
      applyBoost=reconstructISJets(p_cm[0],ShowerHardJets,intrinsic,
				   boostRest,boostNewF);
    if(boostRest.mag()>1.||boostNewF.mag()>1.) return false;
    // final-state reconstruction
    // check if in CMF frame
    Vector3 beta_cm = p_cm[1].findBoostToCM();
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
  catch(Veto) {
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
  if(momConsEq( 0.0, root_s, jets )>0.0) return -1.0;
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
    
    if ( momConsEq( k1, root_s, jets ) < 0.0 ) {
      while ( momConsEq( k2, root_s, jets ) < 0.0 ) {
	k1 = k2; 
	k2 *= 2;       
      }
      while ( fabs( (k1 - k2)/(k1 + k2) ) > 1.e-10 ) { 
	if( momConsEq( k2, root_s, jets ) == 0. ) {
	  return k2; 
	} else {
	  k = (k1+k2)/2.;
	  if ( momConsEq( k, root_s, jets ) > 0 ) {
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
      reconstructTimeLikeJet(child,0);
      // calculate the momentum of the particle
      Lorentz5Momentum pnew=p->momentum()-child->momentum();
      pnew.rescaleMass();
      p->children()[0]->set5Momentum(pnew);
    }
  }
  return emitted;
}

Vector3 QTildeReconstructor::
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
  Vector3 beta = -betam*(k/kp)*oldp.vect();
  // note that (k/kp)*oldp.vect() = oldp.vect()/oldp.vect().mag() but cheaper. 

  // leave this out if it's running properly! 
  if ( betam >= 0 ) return beta;
  else              return Vector3(0., 0., 0.); 
}

bool QTildeReconstructor::
reconstructISJets(Lorentz5Momentum pcm,
		  const vector<ShowerProgenitorPtr> & ShowerHardJets,
		  map<tShowerProgenitorPtr,pair<Energy,double> > intrinsic,
		  Vector3 & boostRest,Vector3 & boostNewF) const {
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
    if(ShowerHardJets[ix]->progenitor()->showerKinematics()) {
      pq.push_back(ShowerHardJets[ix]->progenitor()->
		   showerKinematics()->getBasis()[0]);
    }
    else {
      if (!ShowerHardJets[ix]->original()->parents().empty()) {
	Energy etemp = ShowerHardJets[ix]->original()->
	  parents()[0]->momentum().pz();
	pq.push_back(Lorentz5Momentum(0, 0, etemp, abs(etemp)));
      } 
      else 
	throw Exception() << "QTildeReconstructor::reconstructISJets: "
			   << "Warning, bad pq!!!\n"
			   << Exception::eventerror;
    }
  }
  // add the intrinsic pt if needed
  int iloc=-1;
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
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
  // is this textbook method for solving a quadratic
  // numerically stable if 4AC ~= B^2 ? check Numerical Recipes
  double kp = 1.0;
  Energy2 A = a1*b2*S;
  Energy2 B = sqr(MDY) - (a1*b1+a2*b2)*S - sqr(p1p+p2p);
  Energy2 C = a2*b1*S;
  double rad = 1.-4.*A*C/sqr(B);
  if (rad >= 0) {
    kp = B/(2.*A)*(1.+sqrt(rad));
  }
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
  if (pq[0].z() > 0) {
    beta1 = -beta1; 
    beta2 = -beta2;
  }
  tPVector toBoost;
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    if(!ShowerHardJets[ix]->progenitor()->isFinalState())
      toBoost.push_back(ShowerHardJets[ix]->progenitor());
  }
  // before boost
  Hep3Vector betaboost = Vector3(0, 0, beta1);
  tPPtr parent;
  boostChain(toBoost[0], betaboost,parent);
  if(parent->momentum().e()/pq[0].e()>1.||parent->momentum().z()/pq[0].z()>1.) throw Veto();
  betaboost = Vector3(0, 0, beta2);
  boostChain(toBoost[1], betaboost,parent);
  if(parent->momentum().e()/pq[1].e()>1.||parent->momentum().z()/pq[1].z()>1.) throw Veto();
  boostRest = pcm.findBoostToCM();
  boostNewF = (toBoost[0]->momentum() + toBoost[1]->momentum()).boostVector();
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
    // initial-state radiation
    bool radiated[2] = {false,false};
    ShowerProgenitorPtr initial;
    for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
      // only consider initial-state jets
      if(ShowerHardJets[ix]->progenitor()->isFinalState()) continue;
      initial=ShowerHardJets[ix];
    }
    // if initial state radiation do the reconstruction
    if (!initial->progenitor()->children().empty()) radiated[0] = true;
    Lorentz5Momentum pjet;
      Hep3Vector       boosttorest=-initial->progenitor()->momentum().boostVector();
      Lorentz5Momentum nvect;
      ShowerParticlePtr partner;
      // flag for initial-state restruction procedure
      bool initialrad=radiated[0];
      // check if need to boost to rest frame
      bool gottaBoost = (boosttorest.mag() > 1e-12);
      if(radiated[0]) {
	reconstructDecayJet(initial->progenitor());
	// momentum of decaying particle after ISR
	pjet=initial->progenitor()->momentum()
	  -decay->incomingLines().begin()->second->momentum();
	pjet.rescaleMass();
	Lorentz5Momentum ptest(pjet);
	ptest.boost(-initial->progenitor()->momentum().boostVector());
	// get the n reference vector
	nvect= initial->progenitor()->showerKinematics()->getBasis()[1];
	// find the partner
	partner=initial->progenitor()->
	  partners()[initial->progenitor()->showerKinematics()->
		     splittingFn()->interactionType()];
      }
      // check it is a final state particle and find momentum of the
      // rest of the particles
      Lorentz5Momentum pjetA=pjet;
      int ipartner(-1);
      Lorentz5Momentum p_cm,p_cmnew;
      for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
	// only consider final-state jets
	if(!ShowerHardJets[ix]->progenitor()->isFinalState()) continue;
	// did it radiate
	if(!ShowerHardJets[ix]->progenitor()->children().empty()) {
	  radiated[1] = true;
	  ShowerParticlePtr ptemp=ShowerHardJets[ix]->progenitor()->partners()
	    [ShowerHardJets[ix]->progenitor()->showerKinematics()->
	     splittingFn()->interactionType()];
	  if(ptemp&&!partner) { 
	    if(!ptemp->isFinalState()) {
	      partner=ShowerHardJets[ix]->progenitor();
	      nvect=partner->momentum();
	      nvect.boost(boosttorest);
	      nvect = Lorentz5Momentum(0.,0.5*initial->progenitor()->mass()*
				       nvect.vect().unit());
	      nvect.boost(-boosttorest);
	    }
	  }
	}
	// add momentum
	if(ShowerHardJets[ix]->progenitor()!=partner)
	  p_cm+=ShowerHardJets[ix]->progenitor()->momentum();
	else
	  ipartner=ix;
      }
      p_cm.rescaleMass();
      p_cmnew=p_cm;
      if(partner) initialrad=true;
      // perform the initial-state stuff
      if(initialrad) {
	// now boost the relevant vectors to the rest frame
	if(gottaBoost) {
	  pjet.boost(boosttorest);
	  p_cm.boost(boosttorest);
	  p_cmnew.boost(boosttorest);
	  nvect.boost(boosttorest);
	}
	// reconstruct the partner jet
	Lorentz5Momentum ppartner[2];
	ppartner[0]=partner->momentum();
	if(gottaBoost) ppartner[0].boost(boosttorest);
	reconstructTimeLikeJet(partner,0);
	if(gottaBoost) partner->deepBoost(boosttorest);
	ppartner[1]=partner->momentum();
	// calculate the rescaling parameters
	double k1,k2;
	Lorentz5Momentum qt;
	if(!solveDecayKFactor(initial->progenitor()->mass(),
			      nvect,pjet,p_cm,ppartner,k1,k2,qt)) return false;
	// apply the boosts
	// to the colour partner
	Lorentz5Momentum pnew=ppartner[0];
	pnew *=k1;
	pnew-=qt;
	pnew.setMass(ppartner[1].mass());
	pnew.rescaleEnergy();
	LorentzRotation Trafo=solveBoost(1.,ppartner[1],pnew);
	if(gottaBoost) Trafo.boost(-boosttorest);
	partner->deepTransform(Trafo);
	// and the singlet system
	Trafo = solveBoost(k2, p_cm, p_cmnew);
	p_cmnew.transform(Trafo);
      }
      // check if any final-state radiation
      bool atLeastOnce = radiated[1];
      // collection of pointers to initial hard particle and jet momenta
      // for final boosts
      JetKinVect jetKinematics;
      vector<ShowerProgenitorPtr>::const_iterator cit;
      // boost to the rest frame
      LorentzRotation tboost;
      if(gottaBoost) tboost=LorentzRotation(boosttorest);
      bool finalboost(false);
      Hep3Vector singletboost;
      if(initialrad) {
	singletboost = -p_cm.boostVector();
	finalboost   = singletboost.mag()>1e-12;
	if(finalboost) tboost.boost(singletboost);
	singletboost=-p_cmnew.boostVector();
      }
      int ix=-1;
      // reconstruct the final-state jets and perform boosts
      for(cit = ShowerHardJets.begin(); cit != ShowerHardJets.end(); cit++) {
	  ++ix;
	  if(!(*cit)->progenitor()->isFinalState()) continue;
	  if(ipartner==ix) continue;
	  // only those which don't have initial-state partners
	  JetKinStruct tempJetKin;      
	  tempJetKin.parent = (*cit)->progenitor(); 
	  tempJetKin.p = (*cit)->progenitor()->momentum();
	  if(gottaBoost) tempJetKin.p.transform(tboost);
	  atLeastOnce |= reconstructTimeLikeJet(tempJetKin.parent,0);
	  if(gottaBoost) tempJetKin.parent->deepTransform(tboost); 
	  tempJetKin.q = (*cit)->progenitor()->momentum();
	  jetKinematics.push_back(tempJetKin);
	}
      // find the rescaling factor
      Energy roots=p_cm.m();
      double k = 0.0; 
      if(atLeastOnce) {
	if(jetKinematics.size()>1) {
	    k = solveKfactor(roots, jetKinematics);
	    if(k < 0.) return false;
	  }
	else k=1.;
      }
      // perform the rescaling and boosts
      for(JetKinVect::iterator it = jetKinematics.begin();
	  it != jetKinematics.end(); ++it) {
	LorentzRotation Trafo = LorentzRotation(); 
	// boost for rescaling
	if(atLeastOnce) Trafo = solveBoost(k, it->q, it->p);
	// boost back to lab
	if(finalboost) Trafo.boost(-singletboost);
	if(gottaBoost) Trafo.boost(-boosttorest);
	if(atLeastOnce || gottaBoost || finalboost) it->parent->deepTransform(Trafo);
      }
      Lorentz5Momentum ptotal;
      for(unsigned int ix=0;ix<ShowerHardJets.size();++ix)
	{
	  if(ShowerHardJets[ix]->progenitor()->isFinalState()) 
	    ptotal-=ShowerHardJets[ix]->progenitor()->momentum();
	  else
	    ptotal+=ShowerHardJets[ix]->progenitor()->momentum();
	}
      ptotal-=pjetA;
      ptotal.rescaleMass();
    }
  catch(Veto)
    { return false;}
  return true;
}

bool QTildeReconstructor::
reconstructDecayJet( const tShowerParticlePtr p) const {
  if(p->children().empty()) return false;
  tShowerParticlePtr child;
  // if branching reconstruct time-like child
  child = dynamic_ptr_cast<ShowerParticlePtr>(p->children()[1]);
  if(child)
    {
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
solveDecayKFactor(Energy mb,Lorentz5Momentum n, Lorentz5Momentum pjet, 
		  Lorentz5Momentum pother, Lorentz5Momentum ppartner[2], 
		  double & k1, double & k2,Lorentz5Momentum & qt) const {
  long double d1,d2;
  long double pjn(pjet.vect()*n.vect()),pcn(ppartner[0].vect()*n.vect()),
    pan(pother.vect()*n.vect()),nmag(n.vect().mag2());
  Lorentz5Momentum pn=(pjn/nmag)*n;
  qt=pjet-pn;qt.setE(0.);
  long double pt2=qt.vect().mag2();
  long double Ejet = pjet.e();
  long double pcmag=ppartner[0].vect().mag2();
  long double pamag=pother.vect().mag2();
  long double roots,ea,ec,ds;
  static long double eps=1e-8*GeV;
  d1=1.;
  d2=1.;
  unsigned int ix=0;
  do {
    // calculate new value of d1 using N-R
    d2=-(pjn+d1*pcn)/pan;
    ec = sqrt(sqr(d1)*pcmag+pt2+ppartner[1].mass2());
    ea = sqrt(sqr(d2)*pamag+pother.mass2());
    roots=Ejet+ea+ec;
    ds=d1/ec*pcmag-d2*pamag/ea*pcn/pan;
    d1+=(mb-roots)/ds;
    d2=-(pjn+d1*pcn)/pan;
    // check o.K
    ec = sqrt(sqr(d1)*pcmag+pt2+ppartner[1].mass2());
    ea = sqrt(sqr(d2)*pamag+pother.mass2());
    roots=Ejet+ea+ec;
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
  Vector3 boostv=-pin[0].boostVector();
  for(unsigned int ix=0;ix<pout.size();++ix) {
    pout[ix].boost(boostv);
  }
  double lambda=inverseRescaleingFactor(pout,mon,pin[0].mass());
  // now calculate the p reference vectors 
  for(cit=branchings.begin();cit!=branchings.end();++cit){
    if(!(*cit)->_particle->isFinalState()) continue;
    (*cit)->_p=(*cit)->_particle->momentum();
    (*cit)->_p.boost(boostv);
    (*cit)->_p*=lambda;
    (*cit)->_p.setMass((*cit)->_particle->dataPtr()->mass());
    (*cit)->_p.rescaleEnergy();
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
    Hep3Vector boost=((*cit)->_p+branch->_p).findBoostToCM();
    Lorentz5Momentum pcm = branch->_p;
    pcm.boost(boost);
    (*cit)->_n = Lorentz5Momentum(0.,pcm.vect());
    (*cit)->_n.boost( -boost);
  }
  // now compute the new momenta 
  for(cit=branchings.begin();cit!=branchings.end();++cit){
    if(!(*cit)->_particle->isFinalState()) continue;
    Energy2 dot=(*cit)->_p*(*cit)->_n;
    double beta = 0.5*(sqr((*cit)->_particle->mass())-sqr((*cit)->_p.mass()))/dot;
    Lorentz5Momentum qnew=(*cit)->_p+beta*(*cit)->_n;
    qnew.rescaleMass();
    // compute the boost
    LorentzRotation A=LorentzRotation(boostv);
    LorentzRotation R=solveBoost(qnew,A*(*cit)->_particle->momentum())*A;
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
    Energy sum(0.);
    for(unsigned int ix=0;ix<pout.size();++ix) {
      root[ix] = sqrt(pmag[ix]/sqr(lambda)-sqr(mon[ix]));
      sum+=root[ix];
    }
    // if accuracy reached exit
    if(abs(sum/roots-1.)<1e-10) break;
    // use Newton-Raphson to compute new guess for lambda
    Energy numer(0.),denom(0.);
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

LorentzRotation QTildeReconstructor::solveBoost(const Lorentz5Momentum & newq, 
						const Lorentz5Momentum & oldq ) const {
  Energy k  = oldq.vect().mag();
  Energy q  = newq.vect().mag();
  Energy ek = sqrt(sqr(k)+sqr(oldq.mass()));
  Energy eq = sqrt(sqr(q)+sqr(newq.mass())); 
  double betam = -(k*ek-q*eq)/(sqr(k)+sqr(q)+sqr(newq.mass()));
  Vector3 beta = betam/k*oldq.vect();
  Hep3Vector ax = newq.vect().cross( oldq.vect() ); 
  double delta = newq.vect().angle( oldq.vect() );
  LorentzRotation R;
  if ( ax.mag2()/GeV2 > 1e-16 )
    R.rotate( delta, ax ).boost( beta ); 
  else
    R.boost( beta ); 
  return R;
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
      Energy etemp = (*cit)->_beam->momentum().pz();
      pq.push_back(Lorentz5Momentum(0, 0, etemp, abs(etemp)));
    }
  }
  bool order = (*hard->incoming().begin())->_beam->momentum().pz()/pq[0].z()<0.;
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
  if(pq[0].z()<0) swap(x[0],x[1]);
  double k1=alpha[0]/x[0],k2=beta[1]/x[1];
  double alphanew[2]={alpha[0]/k1,alpha[1]*k2};
  double betanew [2]={beta [0]*k1,beta [1]/k2};
  double boost[2];
  for(unsigned int ix=0;ix<2;++ix) {
    boost[ix] = getBeta(alpha   [ix]+beta   [ix], alpha[ix]   -beta   [ix], 
		       alphanew[ix]+betanew[ix], alphanew[ix]-betanew[ix]);
    if (pq[0].z() > 0) beta[ix]*=-1.;
  }
  // now apply the boosts
  Hep3Vector betaboost = Vector3(0.,0.,boost[0]);
  // apply the boost the the particles
  // first incoming particle
  if(order) swap(pq[0],pq[1]);
  LorentzRotation R;
  R.boost(betaboost);
  branchings=hard->incoming();
  cit=branchings.begin();
  (*cit)->_p=pq[0];
  (*cit)->_n=pq[1];
  (*cit)->setMomenta(R,1.,Lorentz5Momentum());
  // second incoming particle
  betaboost = Vector3(0.,0.,boost[1]);
  R=LorentzRotation(betaboost);
  ++cit;
  (*cit)->_p=pq[1];
  (*cit)->_n=pq[0];
  (*cit)->setMomenta(R,1.,Lorentz5Momentum());
  return true;
}
