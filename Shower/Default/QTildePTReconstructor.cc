// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RelativePtReconstructor class.
//

#include "RelativePtReconstructor.h"
#include "ThePEG/Utilities/Timer.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "RelativePtReconstructor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Kinematics/QtildaShowerKinematics1to2.h"
using namespace Herwig;

void RelativePtReconstructor::persistentOutput(PersistentOStream & os) const {
  os << _massopt;
}

void RelativePtReconstructor::persistentInput(PersistentIStream & is, int) {
  is >> _massopt;
}

ClassDescription<RelativePtReconstructor> RelativePtReconstructor::initRelativePtReconstructor;
// Definition of the static class description member.

void RelativePtReconstructor::Init() {

  static ClassDocumentation<RelativePtReconstructor> documentation
    ("There is no documentation for the RelativePtReconstructor class");


  static Switch<RelativePtReconstructor,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the translation between qtilde and mass",
     &RelativePtReconstructor::_massopt, 0, false, false);
  static SwitchOption interfaceMassOptionMass
    (interfaceMassOption,
     "Mass",
     "Use the basic definition in terms of qtilde",
     0);
  static SwitchOption interfaceMassOptionPt
    (interfaceMassOption,
     "Pt",
     "Use the definition of pt in terms of qtilde and then calculate the mass using pt.",
     1);

}

bool RelativePtReconstructor::
reconstructTimeLikeJet(const tShowerParticlePtr particleJetParent,
		       unsigned int iopt) const {
  if(!particleJetParent)
    {throw Exception() << "must have a particle in Kinematics"
 		       << "Reconstructor::reconstructTimeLikeJet"
 		       << Exception::eventerror;}
  bool emitted=true;
  // if this is not a fixed point in the reconstruction
  if( !(particleJetParent->isReconstructionFixedPoint()) ) 
    {
      // if not a reconstruction fixpoint, dig deeper for all children:
      for ( ParticleVector::const_iterator cit = particleJetParent->children().begin();
 	    cit != particleJetParent->children().end(); ++cit )
	{
	  reconstructTimeLikeJet(dynamic_ptr_cast<ShowerParticlePtr>(*cit),iopt);
	}
    }
  // it is a reconstruction fixpoint, ie kinematical data has to be available 
  else 
    {
      // check if the parent was part of the shower
      ShowerParticlePtr jetGrandParent = dynamic_ptr_cast<ShowerParticlePtr>
 	(particleJetParent->parents()[0]);
      // update if so
      if (jetGrandParent) 
 	{if (jetGrandParent->showerKinematics())
 	    jetGrandParent->showerKinematics()->updateLast(particleJetParent,iopt);}
      // otherwise
      else 
 	{
 	  Energy dm; 
 	  if (particleJetParent->id()==ParticleID::g)
 	    dm=showerVariables()->gluonMass();
 	  else
 	    dm = particleJetParent->data().constituentMass();
 	  if (abs(dm-particleJetParent->momentum().mass())>0.1*MeV
 	      &&particleJetParent->dataPtr()->stable()) 
 	    {
 	      Lorentz5Momentum dum =  particleJetParent->momentum();
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
    {
      particleJetParent->showerKinematics()
 	->updateParent( particleJetParent, particleJetParent->children() );
      Energy pT=particleJetParent->showerKinematics()->pT();
      double zz=particleJetParent->showerKinematics()->z(); 
      Energy mnew=sqrt(sqr(particleJetParent->children()[0]->mass())/zz+
		       sqr(particleJetParent->children()[1]->mass())/(1.-zz)+
		       sqr(pT)/zz/(1.-zz));
      particleJetParent->set5Momentum(Lorentz5Momentum(0.,0.,0.,mnew,mnew));
    }
  return emitted;
}

bool RelativePtReconstructor::
reconstructHardJets(ShowerTreePtr hard) const
{
  //cerr << "testing start of main recon " << endl;
  //cerr << "testing called new B " << endl;
  Timer<1150> timer("KinematicsReconstructor::reconstructHardJets");
  bool radiated[2] = {false,false};
  // find the hard process centre-of-mass energy
  Lorentz5Momentum p_cm[2] = {Lorentz5Momentum(),Lorentz5Momentum()};
  // create a vector of the hard particles
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mit;
  vector<ShowerProgenitorPtr> ShowerHardJets;
  for(mit=hard->incomingLines().begin();mit!=hard->incomingLines().end();++mit)
      ShowerHardJets.push_back((*mit).first);
  for(mit=hard->outgoingLines().begin();mit!=hard->outgoingLines().end();++mit)
      ShowerHardJets.push_back((*mit).first);
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix)
    {
      // final-state jet
      if(ShowerHardJets[ix]->progenitor()->isFinalState())
	{
	  // did it radiate
	  radiated[1] |=ShowerHardJets[ix]->hasEmitted();
	  // add momentum
	  p_cm[1]+=ShowerHardJets[ix]->progenitor()->momentum();
	}
      // initial-state jet
      else
	{
	  // did it radiate
	  radiated[0]|=ShowerHardJets[ix]->hasEmitted();
	  // add momentum
	  p_cm[0]+=ShowerHardJets[ix]->progenitor()->getThePEGBase()->momentum();
	}
    }
  // initial state shuffling
  // the boosts for the initial state
  Vector3 boostRest,boostNewF;
  bool applyBoost(false);
  if(radiated[0])
    applyBoost=reconstructISJets(p_cm[0],ShowerHardJets,boostRest,boostNewF);
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
    if((*cit)->progenitor()->isFinalState()) {
      JetKinStruct tempJetKin;      
      tempJetKin.parent = (*cit)->progenitor(); 
      if(gottaBoost) tempJetKin.parent->boost(beta_cm); 
      tempJetKin.p = (*cit)->progenitor()->momentum();
      atLeastOnce |= reconstructTimeLikeJet((*cit)->progenitor(),0);
      if(!(*cit)->progenitor()->children().empty())
	generateTimeLikeMomenta((*cit)->progenitor(),true);
      tempJetKin.q = (*cit)->progenitor()->momentum();
      jetKinematics.push_back(tempJetKin);  
    }
  }
  // find the rescaling factor
  double k = 0.0; 
  if(atLeastOnce) {
    k = solveKfactor(p_cm[1].mag(), jetKinematics);
    if(k < 0. || k > 1.) return false;
  }
  // perform the rescaling and boosts
  for(JetKinVect::iterator it = jetKinematics.begin();
      it != jetKinematics.end(); ++it) {
    LorentzRotation Trafo = LorentzRotation(); 
    if(atLeastOnce) Trafo = solveBoost(k, it->q, it->p);
    if(gottaBoost) Trafo.boost(-beta_cm);
    if(atLeastOnce || gottaBoost) it->parent->deepTransform(Trafo);
    if(applyBoost)
      {
 	it->parent->deepBoost(boostRest);
 	it->parent->deepBoost(boostNewF);
      }
  }
  return true;
}
void RelativePtReconstructor::
generateTimeLikeMomenta(const tShowerParticlePtr particleJetParent,bool first) const
{

  QtildaShowerKinematics1to2Ptr kin=
    dynamic_ptr_cast<QtildaShowerKinematics1to2Ptr>
    (particleJetParent->showerKinematics());
  Energy pt(kin->pT());
  double phi(kin->phi());
  // generate the momenta as before for the first branching
  if(first)
    {
      Lorentz5Momentum psum;
      for(unsigned int ix=0;ix<particleJetParent->children().size();++ix)
	{
	  ShowerParticlePtr theLast=
	    dynamic_ptr_cast<ShowerParticlePtr>(particleJetParent->children()[ix]);
	  theLast->sudBeta( (sqr(theLast->mass()) + sqr(kin->pT())
			     - sqr( theLast->sudAlpha() )*kin->pVector().m2())
			    / ( 2.*theLast->sudAlpha()*kin->p_dot_n() ) );
	  Energy px(pt*cos(phi)),py(pt*sin(phi));
	  if(ix!=0){px*=-1.;py*=-1.;}
	  theLast->set5Momentum(kin->sudakov2Momentum(theLast->sudAlpha(),
						      theLast->sudBeta(), 
						      px,py,0));
	  psum+=theLast->momentum();
	  if(!theLast->children().empty())
	    generateTimeLikeMomenta(theLast,false);
	}
      psum.setMass(particleJetParent->mass());
      particleJetParent->set5Momentum(psum);
    }
  else
    {
      // compute the p reference vector
      Energy A(particleJetParent->momentum().e()+
	       particleJetParent->momentum().vect().mag()),
	m(particleJetParent->dataPtr()->mass());
      Energy modp(0.5*(A+m)*(A-m)/A);
      Lorentz5Momentum p(m,modp*particleJetParent->momentum().vect().unit());
      // n reference vector
      Lorentz5Momentum n(0.,-particleJetParent->momentum().vect().unit());
      Lorentz5Momentum ptest(p+(particleJetParent->mass()-m)*
			     (particleJetParent->mass()+m)/2./(p*n)*n);
      ptest.rescaleMass();
      double alphaold(particleJetParent->sudAlpha());
      Energy pdotn=p*n;
      Lorentz5Momentum psum;
      for(unsigned int ix=0;ix<particleJetParent->children().size();++ix)
	{
	  ShowerParticlePtr theLast=
	    dynamic_ptr_cast<ShowerParticlePtr>(particleJetParent->children()[ix]);
	  double alpha(theLast->sudAlpha()/alphaold);
	  double beta=0.5*(sqr(theLast->mass())+sqr(pt)-sqr(alpha*m))/alpha/pdotn;
	  Energy px(pt*cos(phi)),py(pt*sin(phi));
	  if(ix!=0){px*=-1.;py*=-1.;}
	  Lorentz5Momentum pnew(px,py,alpha*modp-beta,
				alpha*sqrt(sqr(modp)+sqr(m))+beta);
	  pnew.rotateUz(p.vect().unit());
	  pnew.rescaleMass();
	  theLast->set5Momentum(pnew);
	  psum+=pnew;
	  if(!theLast->children().empty())
	    generateTimeLikeMomenta(theLast,false);
	}
    }
}
