// -*- C++ -*-
//
// FS_QTildeShowerKinematics1to2.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_QTildeShowerKinematics1to2 class.
//

#include "FS_QTildeShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"

using namespace Herwig;

void FS_QTildeShowerKinematics1to2::
updateParameters(tShowerParticlePtr theParent,
		 tShowerParticlePtr theChild0,
		 tShowerParticlePtr theChild1,
		 bool setAlpha) const {
  const ShowerParticle::Parameters & parent = theParent->showerParameters();
  ShowerParticle::Parameters & child0 = theChild0->showerParameters();
  ShowerParticle::Parameters & child1 = theChild1->showerParameters();
  // determine alphas of children according to interpretation of z
  if ( setAlpha ) {
    child0.alpha =      z() * parent.alpha;
    child1.alpha = (1.-z()) * parent.alpha;
  }
  // set the values
  double cphi = cos(phi());
  double sphi = sin(phi());

  child0.ptx =  pT() * cphi +     z() * parent.ptx;
  child0.pty =  pT() * sphi +     z() * parent.pty;
  child0.pt  = sqrt( sqr(child0.ptx) + sqr(child0.pty) );

  child1.ptx = -pT() * cphi + (1.-z())* parent.ptx;
  child1.pty = -pT() * sphi + (1.-z())* parent.pty;
  child1.pt  = sqrt( sqr(child1.ptx) + sqr(child1.pty) );
}



void FS_QTildeShowerKinematics1to2::
updateChildren(const tShowerParticlePtr theParent, 
	       const ShowerParticleVector & theChildren,
	       bool angularOrder) const {
  if(theChildren.size() != 2)
    throw Exception() << "FS_QTildeShowerKinematics1to2::updateChildren() "
		      << "Warning! too many children!" << Exception::eventerror;
  // copy scales etc
  Energy dqtilde = scale();
  // note that 1st child gets z, 2nd gets (1-z) by our convention.
  if(angularOrder) {
    theChildren[0]->setEvolutionScale(    z() * dqtilde);
    theChildren[1]->setEvolutionScale((1.-z())* dqtilde);
  }
  else {
    theChildren[0]->setEvolutionScale(dqtilde);
    theChildren[1]->setEvolutionScale(dqtilde);
  }
  updateParameters(theParent, theChildren[0], theChildren[1], true);
  // set up the colour connections
  splittingFn()->colourConnection(theParent,theChildren[0],theChildren[1],false);
  // make the products children of the parent
  theParent->addChild(theChildren[0]);
  theParent->addChild(theChildren[1]);
}

void FS_QTildeShowerKinematics1to2::
reconstructParent(const tShowerParticlePtr theParent, 
	     const ParticleVector & theChildren ) const {
  if(theChildren.size() != 2) 
    throw Exception() << "FS_QTildeShowerKinematics1to2::updateParent() " 
		      << "Warning! too many children!" 
		      << Exception::eventerror;
  ShowerParticlePtr c1 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[0]);
  ShowerParticlePtr c2 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[1]);
  theParent->showerParameters().beta= 
    c1->showerParameters().beta + c2->showerParameters().beta; 
  theParent->set5Momentum( c1->momentum() + c2->momentum() );
}

void FS_QTildeShowerKinematics1to2::reconstructLast(const tShowerParticlePtr theLast,
						    unsigned int iopt, 
						    Energy mass) const {
  // set beta component and consequently all missing data from that,
  // using the nominal (i.e. PDT) mass.
  Energy theMass = mass > ZERO  ?  mass : theLast->data().constituentMass();
  ShowerParticle::Parameters & last = theLast->showerParameters();
  last.beta = ( sqr(theMass) + sqr(last.pt) - sqr(last.alpha) * pVector().m2() )
    / ( 2. * last.alpha * p_dot_n() );
  // set that new momentum
  theLast->set5Momentum(sudakov2Momentum( last.alpha, last.beta, 
					  last.ptx, last.pty, iopt) );
}

void FS_QTildeShowerKinematics1to2::initialize(ShowerParticle & particle,PPtr) {
  // set the basis vectors
  Lorentz5Momentum p,n;
  if(particle.perturbative()!=0) {
    // find the partner and its momentum
    ShowerParticlePtr partner=particle.partner();
    Lorentz5Momentum ppartner(partner->momentum());
    // momentum of the emitting particle
    p = particle.momentum();
    Lorentz5Momentum pcm;
    // if the partner is a final-state particle then the reference
    // vector is along the partner in the rest frame of the pair
    if(partner->isFinalState()) {
      Boost boost=(p + ppartner).findBoostToCM();
      pcm = ppartner;
      pcm.boost(boost);
      n = Lorentz5Momentum(ZERO,pcm.vect());
      n.boost( -boost);
    }
    else if(!partner->isFinalState()) {
      // if the partner is an initial-state particle then the reference
      // vector is along the partner which should be massless
      if(particle.perturbative()==1)
	{n = Lorentz5Momentum(ZERO,ppartner.vect());}
      // if the partner is an initial-state decaying particle then the reference
      // vector is along the backwards direction in rest frame of decaying particle
      else {
	Boost boost=ppartner.findBoostToCM();
	pcm = p;
	pcm.boost(boost);
	n = Lorentz5Momentum( ZERO, -pcm.vect()); 
	n.boost( -boost);
      } 
    } 
  }
  else if(particle.initiatesTLS()) {
    tShoKinPtr kin=dynamic_ptr_cast<ShowerParticlePtr>
      (particle.parents()[0]->children()[0])->showerKinematics();
    p = kin->getBasis()[0];
    n = kin->getBasis()[1];
  }
  else  {
    tShoKinPtr kin=dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
      ->showerKinematics();
    p = kin->getBasis()[0];
    n = kin->getBasis()[1];
  }
  // set the basis vectors
  setBasis(p,n);
}

void FS_QTildeShowerKinematics1to2::updateParent(const tShowerParticlePtr parent, 
						 const ShowerParticleVector & children,
						 bool) const {
  IdList ids(3);
  ids[0] = parent->id();
  ids[1] = children[0]->id();
  ids[2] = children[1]->id();
  const vector<Energy> & virtualMasses = SudakovFormFactor()->virtualMasses(ids);
  if(children[0]->children().empty()) children[0]->setVirtualMass(virtualMasses[1]);
  if(children[1]->children().empty()) children[1]->setVirtualMass(virtualMasses[2]);
  // compute the new pT of the branching
  Energy2 pt2=sqr(z()*(1.-z()))*sqr(scale())
    - sqr(children[0]->virtualMass())*(1.-z())
    - sqr(children[1]->virtualMass())*    z() ;
  if(ids[0]!=ParticleID::g) pt2 += z()*(1.-z())*sqr(virtualMasses[0]);
  Energy2 q2 = 
    sqr(children[0]->virtualMass())/z() + 
    sqr(children[1]->virtualMass())/(1.-z()) +
    pt2/z()/(1.-z());
  if(pt2<ZERO) {
    parent->setVirtualMass(ZERO);
  }
  else {
    parent->setVirtualMass(sqrt(q2));
    pT(sqrt(pt2));
  }
}

void FS_QTildeShowerKinematics1to2::
resetChildren(const tShowerParticlePtr theParent, 
	      const ShowerParticleVector & theChildren) const {
  updateParameters(theParent, theChildren[0], theChildren[1], false);

  for(unsigned int ix=0;ix<theChildren.size();++ix) {
    if(theChildren[ix]->children().empty()) continue;
    ShowerParticleVector children;
    for(unsigned int iy=0;iy<theChildren[ix]->children().size();++iy)
      children.push_back(dynamic_ptr_cast<ShowerParticlePtr>
			 (theChildren[ix]->children()[iy]));
    theChildren[ix]->showerKinematics()->resetChildren(theChildren[ix],children);
  }
}
