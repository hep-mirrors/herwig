// -*- C++ -*-
//
// Decay_QTildeShowerKinematics1to2.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Decay_QTildeShowerKinematics1to2 class.
//

#include "Decay_QTildeShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

void Decay_QTildeShowerKinematics1to2::
updateChildren(const tShowerParticlePtr theParent, 
	       const ShowerParticleVector & theChildren,
	       bool angularOrder ) const {
  if(theChildren.size() != 2)
    throw Exception() <<  "Decay_QTildeShowerKinematics1to2::updateChildren() " 
 		      << "Warning! too many children!" << Exception::eventerror;
  // copy scales etc
  Energy dqtilde = scale();
  double dz = z(); 
  double dphi = phi();
  // set the values
  // if(theParent->showerParameters().alpha==0.0) {
  //   theParent->showerParameters().alpha=1.;
  // }
  if(angularOrder) {
    theChildren[0]->setEvolutionScale(        dqtilde);
    theChildren[1]->setEvolutionScale((1.-dz)*dqtilde);
  }
  else {
    theChildren[0]->setEvolutionScale(        dqtilde);
    theChildren[1]->setEvolutionScale(        dqtilde);
  }
  // determine alphas of children according to interpretation of z
  theChildren[0]->showerParameters().alpha=    dz *theParent->showerParameters().alpha; 
  theChildren[1]->showerParameters().alpha=(1.-dz)*theParent->showerParameters().alpha;
  theChildren[0]->showerParameters().px=   pT()*cos(dphi) +     
    dz *theParent->showerParameters().px;
  theChildren[0]->showerParameters().py=   pT()*sin(dphi) + 
    dz *theParent->showerParameters().py;
  theChildren[1]->showerParameters().px= - pT()*cos(dphi) + 
    (1.-dz)*theParent->showerParameters().px;
  theChildren[1]->showerParameters().py= - pT()*sin(dphi) + 
    (1.-dz)*theParent->showerParameters().py;
  for(unsigned int ix=0;ix<2;++ix)
    theChildren[ix]->showerParameters().pt=
      sqrt(sqr(theChildren[ix]->showerParameters().px)+
	   sqr(theChildren[ix]->showerParameters().py));
  // set up the colour connections
  splittingFn()->colourConnection(theParent,theChildren[0],theChildren[1],false);
  // make the products children of the parent
  theParent->addChild(theChildren[0]);
  theParent->addChild(theChildren[1]);
}

void Decay_QTildeShowerKinematics1to2::
reconstructParent( const tShowerParticlePtr, const ParticleVector &) const {
  throw Exception() << "Decay_QTildeShowerKinematics1to2::updateParent not implemented"
		    << Exception::abortnow;
}

void Decay_QTildeShowerKinematics1to2::
reconstructLast(const tShowerParticlePtr theLast,
		unsigned int iopt,Energy mass) const {
  // set beta component and consequently all missing data from that,
  // using the nominal (i.e. PDT) mass.
  Energy theMass = mass>ZERO ? mass : theLast->data().constituentMass(); 
  theLast->showerParameters().beta=
    (sqr(theMass) + sqr(theLast->showerParameters().pt)
     - sqr( theLast->showerParameters().alpha )*pVector().m2())
    / ( 2.*theLast->showerParameters().alpha*p_dot_n() );   
  // set that new momentum  
  theLast->set5Momentum(  sudakov2Momentum( theLast->showerParameters().alpha, 
					    theLast->showerParameters().beta, 
					    theLast->showerParameters().px,
					    theLast->showerParameters().py,iopt));
}

void Decay_QTildeShowerKinematics1to2::initialize(ShowerParticle & particle,PPtr) {
  Lorentz5Momentum p, n, ppartner, pcm;
  assert(particle.perturbative()!=1);
  // this is for the initial decaying particle
  if(particle.perturbative()==2) {
    p = particle.momentum();
    ShowerParticlePtr partner=particle.partner();
    Lorentz5Momentum ppartner(partner->momentum());
    // reomved to make inverse recon work properly
    //if(partner->getThePEGBase()) ppartner=partner->getThePEGBase()->momentum();
    pcm=ppartner;
    Boost boost(p.findBoostToCM());
    pcm.boost(boost);
    n = Lorentz5Momentum( ZERO,0.5*p.mass()*pcm.vect().unit()); 
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
