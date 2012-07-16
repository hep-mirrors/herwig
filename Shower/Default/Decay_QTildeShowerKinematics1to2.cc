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
updateChildren(const tShowerParticlePtr parent, 
	       const ShowerParticleVector & children,
	       ShowerPartnerType::Type partnerType) const {
  assert(children.size() == 2);
  // set the values
  if(parent->showerVariables().empty()) {
    parent->showerVariables().resize(3);
    parent->showerParameters().resize(2);
    parent->showerParameters()[0]=1.;
  }
  for(unsigned int ix=0;ix<2;++ix) {
    children[ix]->showerVariables() .resize(3);
    children[ix]->showerParameters().resize(2);
  }
  // calculate the scales
  splittingFn()->evaluateDecayScales(partnerType,scale(),z(),parent,
				     children[0],children[1]);
  // determine alphas of children according to interpretation of z
  children[0]->showerParameters()[0]=    z() *parent->showerParameters()[0]; 
  children[1]->showerParameters()[0]=(1.-z())*parent->showerParameters()[0];
  children[0]->showerVariables()[0]=   pT()*cos(phi()) +     
    z() *parent->showerVariables()[0];
  children[0]->showerVariables()[1]=   pT()*sin(phi()) + 
    z() *parent->showerVariables()[1];
  children[1]->showerVariables()[0]= - pT()*cos(phi()) + 
    (1.-z())*parent->showerVariables()[0];
  children[1]->showerVariables()[1]= - pT()*sin(phi()) + 
    (1.-z())*parent->showerVariables()[1];
  for(unsigned int ix=0;ix<2;++ix)
    children[ix]->showerVariables()[2]=
      sqrt(sqr(children[ix]->showerVariables()[0])+
	   sqr(children[ix]->showerVariables()[1]));
  // set up the colour connections
  splittingFn()->colourConnection(parent,children[0],children[1],partnerType,false);
  // make the products children of the parent
  parent->addChild(children[0]);
  parent->addChild(children[1]);
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

void Decay_QTildeShowerKinematics1to2::initialize(ShowerParticle & particle,PPtr) {
  Lorentz5Momentum p, n, ppartner, pcm;
  assert(particle.perturbative()!=1);
  // this is for the initial decaying particle
  if(particle.perturbative()==2) {
    p = particle.momentum();
    ShowerParticlePtr partner=particle.partner();
    Lorentz5Momentum ppartner(partner->momentum());
    // reomved to make inverse recon work properly
    //if(partner->thePEGBase()) ppartner=partner->thePEGBase()->momentum();
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
