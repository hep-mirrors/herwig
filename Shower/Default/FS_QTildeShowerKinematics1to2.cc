// -*- C++ -*-
//
// FS_QTildeShowerKinematics1to2.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_QTildeShowerKinematics1to2 class.
//

#include "FS_QTildeShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Shower/SplittingFunctions/SplittingFunction.h"
#include "Herwig/Shower/Base/ShowerParticle.h"
#include "ThePEG/Utilities/Debug.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Shower/Base/Evolver.h"
#include "Herwig/Shower/Base/PartnerFinder.h"
#include "Herwig/Shower/Base/ShowerModel.h"
#include "Herwig/Shower/Base/KinematicsReconstructor.h"
#include "Herwig/Shower/Base/ShowerVertex.h"

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
updateChildren(const tShowerParticlePtr parent, 
	       const ShowerParticleVector & children,
	       ShowerPartnerType::Type partnerType) const {
  assert(children.size()==2);
  // calculate the scales
  splittingFn()->evaluateFinalStateScales(partnerType,scale(),z(),parent,
					  children[0],children[1]);
  // update the parameters
  updateParameters(parent, children[0], children[1], true);
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
  // sort out the helicity stuff 
  if(! ShowerHandler::currentHandler()->evolver()->correlations()) return;
  SpinPtr pspin(parent->spinInfo());
  if(!pspin ||  !ShowerHandler::currentHandler()->evolver()->spinCorrelations() ) return;
  Energy2 t = sqr(scale())*z()*(1.-z());
  IdList ids;
  ids.push_back(parent->id());
  ids.push_back(children[0]->id());
  ids.push_back(children[1]->id());
  // create the vertex
  SVertexPtr vertex(new_ptr(ShowerVertex()));
  // set the matrix element
  vertex->ME(splittingFn()->matrixElement(z(),t,ids,phi(),true));
  // set the incoming particle for the vertex
  parent->spinInfo()->decayVertex(vertex);
  for(ShowerParticleVector::const_iterator pit=children.begin();
      pit!=children.end();++pit) {
    // construct the spin info for the children
    constructSpinInfo(*pit,true);
    // connect the spinInfo object to the vertex
    (*pit)->spinInfo()->productionVertex(vertex);
  }
}

void FS_QTildeShowerKinematics1to2::
reconstructParent(const tShowerParticlePtr parent, 
		  const ParticleVector & children ) const {
  assert(children.size() == 2);
  ShowerParticlePtr c1 = dynamic_ptr_cast<ShowerParticlePtr>(children[0]);
  ShowerParticlePtr c2 = dynamic_ptr_cast<ShowerParticlePtr>(children[1]);
  parent->showerParameters().beta= 
    c1->showerParameters().beta + c2->showerParameters().beta; 
  Lorentz5Momentum pnew = c1->momentum() + c2->momentum();
  Energy2 m2 = sqr(pT())/z()/(1.-z()) + sqr(c1->mass())/z()
    + sqr(c2->mass())/(1.-z());
  pnew.setMass(sqrt(m2));
  parent->set5Momentum( pnew );
}

void FS_QTildeShowerKinematics1to2::reconstructLast(const tShowerParticlePtr last,
						    Energy mass) const {
  // set beta component and consequently all missing data from that,
  // using the nominal (i.e. PDT) mass.
  Energy theMass = mass > ZERO  ?  mass : last->data().constituentMass();
  ShowerParticle::Parameters & lastParam = last->showerParameters();
  Energy2 denom = 2. * lastParam.alpha * p_dot_n();
  if(abs(denom)/(sqr(pVector().e())+pVector().rho2())<1e-10) {
    throw KinematicsReconstructionVeto();
  }
  lastParam.beta = ( sqr(theMass) + sqr(lastParam.pt) - sqr(lastParam.alpha) * pVector().m2() )
    / denom;
  // set that new momentum
  Lorentz5Momentum newMomentum = sudakov2Momentum( lastParam.alpha, lastParam.beta,
						   lastParam.ptx  , lastParam.pty);
  newMomentum.setMass(theMass);
  newMomentum.rescaleEnergy();
  if(last->data().stable()) {
    last->set5Momentum( newMomentum );
  }
  else {
    last->boost(last->momentum().findBoostToCM());
    last->boost(newMomentum.boostVector());
  }
}

void FS_QTildeShowerKinematics1to2::initialize(ShowerParticle & particle,PPtr) {
  // set the basis vectors
  Lorentz5Momentum p,n;
  Frame frame;
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
    frame = BackToBack;
  }
  else if(particle.initiatesTLS()) {
    tShoKinPtr kin=dynamic_ptr_cast<ShowerParticlePtr>
      (particle.parents()[0]->children()[0])->showerKinematics();
    p = kin->getBasis()[0];
    n = kin->getBasis()[1];
    frame = kin->frame();
  }
  else  {
    tShoKinPtr kin=dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
      ->showerKinematics();
    p = kin->getBasis()[0];
    n = kin->getBasis()[1];
    frame = kin->frame();
  }
  // set the basis vectors
  setBasis(p,n,frame);
}

void FS_QTildeShowerKinematics1to2::updateParent(const tShowerParticlePtr parent, 
						 const ShowerParticleVector & children,
						 ShowerPartnerType::Type) const {
  IdList ids(3);
  ids[0] = parent->id();
  ids[1] = children[0]->id();
  ids[2] = children[1]->id();
  const vector<Energy> & virtualMasses = SudakovFormFactor()->virtualMasses(ids);
  if(children[0]->children().empty()) children[0]->virtualMass(virtualMasses[1]);
  if(children[1]->children().empty()) children[1]->virtualMass(virtualMasses[2]);
  // compute the new pT of the branching
  Energy2 pt2=sqr(z()*(1.-z()))*sqr(scale())
    - sqr(children[0]->virtualMass())*(1.-z())
    - sqr(children[1]->virtualMass())*    z() ;
  if(ids[0]!=ParticleID::g) pt2 += z()*(1.-z())*sqr(virtualMasses[0]);
  if(pt2>ZERO) {
    Energy2 q2 = 
      sqr(children[0]->virtualMass())/z() + 
      sqr(children[1]->virtualMass())/(1.-z()) +
      pt2/z()/(1.-z());
    parent->virtualMass(sqrt(q2));
    pT(sqrt(pt2));
  }
  else {
    parent->virtualMass(ZERO);
  }
}

void FS_QTildeShowerKinematics1to2::
resetChildren(const tShowerParticlePtr parent, 
	      const ShowerParticleVector & children) const {
  updateParameters(parent, children[0], children[1], false);
  for(unsigned int ix=0;ix<children.size();++ix) {
    if(children[ix]->children().empty()) continue;
    ShowerParticleVector newChildren;
    for(unsigned int iy=0;iy<children[ix]->children().size();++iy)
      newChildren.push_back(dynamic_ptr_cast<ShowerParticlePtr>
			    (children[ix]->children()[iy]));
    children[ix]->showerKinematics()->resetChildren(children[ix],newChildren);
  }
}
