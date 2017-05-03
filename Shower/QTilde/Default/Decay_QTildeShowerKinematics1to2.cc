// -*- C++ -*-
//
// Decay_QTildeShowerKinematics1to2.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Decay_QTildeShowerKinematics1to2 class.
//

#include "Decay_QTildeShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Shower/Core/SplittingFunctions/SplittingFunction.h"
#include "Herwig/Shower/Core/Base/ShowerParticle.h"
#include <cassert>
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Shower/Core/Base/ShowerVertex.h"

using namespace Herwig;

void Decay_QTildeShowerKinematics1to2::
updateChildren(const tShowerParticlePtr parent, 
	       const ShowerParticleVector & children,
	       ShowerPartnerType partnerType,
	       bool massVeto) const {
  assert(children.size() == 2);
  // calculate the scales
  splittingFn()->evaluateDecayScales(partnerType,scale(),z(),parent,
				     children[0],children[1]);
  // set the maximum virtual masses
  IdList ids(3);
  ids[0] = parent->dataPtr();
  ids[1] = children[0]->dataPtr();
  ids[2] = children[1]->dataPtr();
  const vector<Energy> & virtualMasses = SudakovFormFactor()->virtualMasses(ids);
  Energy2 q2 = sqr(virtualMasses[0])-(1.-z())*sqr(scale());
  children[0]->virtualMass(sqrt(q2));
  if(massVeto) {
    children[1]->scales().Max_Q2 = (1.-z())/z()*(z()*sqr(virtualMasses[0])-q2);
  }
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
    (**pit).showerBasis(parent->showerBasis(),true);
    (**pit).setShowerMomentum(true);
  }
}

void Decay_QTildeShowerKinematics1to2::
reconstructParent( const tShowerParticlePtr, const ParticleVector &) const {
  throw Exception() << "Decay_QTildeShowerKinematics1to2::reconstructParent not implemented"
		    << Exception::abortnow;
}

void Decay_QTildeShowerKinematics1to2::
reconstructLast(const tShowerParticlePtr last, Energy mass) const {
  // set beta component and consequently all missing data from that,
  // using the nominal (i.e. PDT) mass.
  Energy theMass = mass > ZERO ? mass : last->data().constituentMass(); 
  last->showerParameters().beta=
    (sqr(theMass) + sqr(last->showerParameters().pt)
     - sqr( last->showerParameters().alpha )*last->showerBasis()->pVector().m2())
    / ( 2.*last->showerParameters().alpha*last->showerBasis()->p_dot_n() );
  // set that new momentum  
  last->set5Momentum( last->showerBasis()->sudakov2Momentum( last->showerParameters().alpha, 
							     last->showerParameters().beta, 
							     last->showerParameters().ptx,
							     last->showerParameters().pty) );
}

void Decay_QTildeShowerKinematics1to2::updateParent(const tShowerParticlePtr parent, 
						    const ShowerParticleVector & children,
						    ShowerPartnerType) const {
  IdList ids(3);
  ids[0] = parent->dataPtr();
  ids[1] = children[0]->dataPtr();
  ids[2] = children[1]->dataPtr();
  const vector<Energy> & virtualMasses = SudakovFormFactor()->virtualMasses(ids);
  children[0]->virtualMass(sqrt(sqr(virtualMasses[0])-(1.-z())*sqr(scale())));
  if(children[1]->children().empty()) children[1]->virtualMass(virtualMasses[2]);
  // compute the new pT of the branching
  Energy2 pt2=(1.-z())*(z()*sqr(virtualMasses[0])-sqr(children[0]->virtualMass()))
    -z()*sqr(children[1]->virtualMass());
  if(pt2>ZERO) {
    pT(sqrt(pt2));
  }
  else {
    parent->virtualMass(ZERO);
  }
}
