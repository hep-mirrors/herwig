// -*- C++ -*-
//
// Decay_QTildeShowerKinematics1to2.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Decay_QTildeShowerKinematics1to2 class.
//

#include "Decay_QTildeShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Shower/SplittingFunctions/SplittingFunction.h"
#include "Herwig/Shower/Base/ShowerParticle.h"
#include <cassert>
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Shower/Base/Evolver.h"
#include "Herwig/Shower/Base/PartnerFinder.h"
#include "Herwig/Shower/Base/ShowerModel.h"
#include "Herwig/Shower/Base/KinematicsReconstructor.h"
#include "Herwig/Shower/Base/ShowerVertex.h"

using namespace Herwig;

void Decay_QTildeShowerKinematics1to2::
updateChildren(const tShowerParticlePtr parent, 
	       const ShowerParticleVector & children,
	       ShowerPartnerType::Type partnerType) const {
  assert(children.size() == 2);
  // calculate the scales
  splittingFn()->evaluateDecayScales(partnerType,scale(),z(),parent,
				     children[0],children[1]);
  // determine alphas of children according to interpretation of z
  const ShowerParticle::Parameters & params = parent->showerParameters();
  ShowerParticle::Parameters & child0 = children[0]->showerParameters();
  ShowerParticle::Parameters & child1 = children[1]->showerParameters();
  child0.alpha =     z()  * params.alpha; 
  child1.alpha = (1.-z()) * params.alpha;

  child0.ptx =  pT() * cos(phi()) +     z()* params.ptx;
  child0.pty =  pT() * sin(phi()) +     z()* params.pty;
  child0.pt  = sqrt( sqr(child0.ptx) + sqr(child0.pty) );

  child1.ptx = -pT() * cos(phi()) + (1.-z()) * params.ptx;
  child1.pty = -pT() * sin(phi()) + (1.-z()) * params.pty;
  child1.pt  = sqrt( sqr(child1.ptx) + sqr(child1.pty) );

  // set up the colour connections
  splittingFn()->colourConnection(parent,children[0],children[1],partnerType,false);
  // make the products children of the parent
  parent->addChild(children[0]);
  parent->addChild(children[1]);
  // set the momenta of the children
  for(ShowerParticleVector::const_iterator pit=children.begin();
      pit!=children.end();++pit) {
    setMomentum(*pit,true);
  }
}

void Decay_QTildeShowerKinematics1to2::
reconstructParent( const tShowerParticlePtr, const ParticleVector &) const {
  throw Exception() << "Decay_QTildeShowerKinematics1to2::updateParent not implemented"
		    << Exception::abortnow;
}

void Decay_QTildeShowerKinematics1to2::
reconstructLast(const tShowerParticlePtr last, Energy mass) const {
  // set beta component and consequently all missing data from that,
  // using the nominal (i.e. PDT) mass.
  Energy theMass = mass>ZERO ? mass : last->data().constituentMass(); 
  last->showerParameters().beta=
    (sqr(theMass) + sqr(last->showerParameters().pt)
     - sqr( last->showerParameters().alpha )*pVector().m2())
    / ( 2.*last->showerParameters().alpha*p_dot_n() );   
  // set that new momentum  
  last->set5Momentum( sudakov2Momentum( last->showerParameters().alpha, 
					last->showerParameters().beta, 
					last->showerParameters().ptx,
					last->showerParameters().pty) );
}

void Decay_QTildeShowerKinematics1to2::initialize(ShowerParticle & particle,PPtr) {
  Lorentz5Momentum p, n, ppartner, pcm;
  Frame frame;
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
    frame = Rest;
  }
  else {
    tShoKinPtr kin=dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
      ->showerKinematics();
    p = kin->getBasis()[0];
    n = kin->getBasis()[1];
    frame = kin->frame();
  }
  setBasis(p,n,frame);
}
