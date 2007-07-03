// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Decay_QtildaShowerKinematics1to2 class.
//

#include "Decay_QtildaShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

void Decay_QtildaShowerKinematics1to2::
updateChildren(const tShowerParticlePtr theParent, 
	       const ShowerParticleVector theChildren ) const {
  if(theChildren.size() != 2)
    throw Exception() <<  "Decay_QtildaShowerKinematics1to2::updateChildren() " 
 		      << "Warning! too many children!" << Exception::eventerror;
  // get the interaction type
  const ShowerIndex::InteractionType interaction =splittingFn()->interactionType();
  // copy scales etc
  Energy dqtilde = scale();
  double dz = z(); 
  double dphi = phi();
  // set the values
  if(theParent->showerVariables().empty()) {
    theParent->showerVariables().resize(3);
    theParent->showerParameters().resize(2);
    theParent->showerParameters()[0]=1.;
  }
  for(unsigned int ix=0;ix<2;++ix) {
    theChildren[ix]->showerVariables() .resize(3);
    theChildren[ix]->showerParameters().resize(2);
  }
  theChildren[0]->setEvolutionScale(interaction,         dqtilde);
  theChildren[1]->setEvolutionScale(interaction, (1.-dz)*dqtilde);
  // determine alphas of children according to interpretation of z
  theChildren[0]->showerParameters()[0]=    dz *theParent->showerParameters()[0]; 
  theChildren[1]->showerParameters()[0]=(1.-dz)*theParent->showerParameters()[0];
  theChildren[0]->showerVariables()[0]=   pT()*cos(dphi) +     
    dz *theParent->showerVariables()[0];
  theChildren[0]->showerVariables()[1]=   pT()*sin(dphi) + 
    dz *theParent->showerVariables()[1];
  theChildren[1]->showerVariables()[0]= - pT()*cos(dphi) + 
    (1.-dz)*theParent->showerVariables()[0];
  theChildren[1]->showerVariables()[1]= - pT()*sin(dphi) + 
    (1.-dz)*theParent->showerVariables()[1];
  for(unsigned int ix=0;ix<2;++ix)
    theChildren[ix]->showerVariables()[2]=
      sqrt(sqr(theChildren[ix]->showerVariables()[0])+
	   sqr(theChildren[ix]->showerVariables()[1]));
  // set up the colour connections
  splittingFn()->colourConnection(theParent,theChildren[0],theChildren[1],false);
  // make the products children of the parent
  theParent->addChild(theChildren[0]);
  theParent->addChild(theChildren[1]);
}

void Decay_QtildaShowerKinematics1to2::
reconstructParent( const tShowerParticlePtr, const ParticleVector ) const {
  throw Exception() << "Decay_QtildaShowerKinematics1to2::updateParent not implemented"
		    << Exception::abortnow;
}

void Decay_QtildaShowerKinematics1to2::reconstructLast(const tShowerParticlePtr theLast,
						  unsigned int iopt) const {
  // set beta component and consequently all missing data from that,
  // using the nominal (i.e. PDT) mass.
  Energy theMass = theLast->data().constituentMass(); 
  theLast->showerParameters()[1]=
    (sqr(theMass) + sqr(theLast->showerVariables()[2])
     - sqr( theLast->showerParameters()[0] )*pVector().m2())
    / ( 2.*theLast->showerParameters()[0]*p_dot_n() );   
  // set that new momentum  
  theLast->set5Momentum(  sudakov2Momentum( theLast->showerParameters()[0], 
					    theLast->showerParameters()[1], 
					    theLast->showerVariables()[0],
					    theLast->showerVariables()[1],iopt));
}

void Decay_QtildaShowerKinematics1to2::initialize(ShowerParticle & particle,PPtr) {
  Lorentz5Momentum p, n, ppartner, pcm;
  assert(particle.perturbative()!=1);
  // this is for the initial decaying particle
  if(particle.perturbative()==2) {
    p = particle.momentum();
    ShowerParticlePtr partner=particle.partners()[splittingFn()->interactionType()];
    Lorentz5Momentum ppartner(partner->momentum());
    if(partner->getThePEGBase()) ppartner=partner->getThePEGBase()->momentum();
    pcm=ppartner;
    Boost boost(p.findBoostToCM());
    pcm.boost(boost);
    n = Lorentz5Momentum( 0.0*MeV,0.5*p.mass()*pcm.vect().unit()); 
    n.boost( -boost);
  }
  else {
    tShoKinPtr kin=dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
      ->showerKinematics();
    p = kin->getBasis()[0];
    n = kin->getBasis()[1];
  }
  setBasis(p,n);
}
