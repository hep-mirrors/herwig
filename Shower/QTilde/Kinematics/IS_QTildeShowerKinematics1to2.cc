// -*- C++ -*-
//
// IS_QTildeShowerKinematics1to2.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_QTildeShowerKinematics1to2 class.
//

#include "IS_QTildeShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"
#include "ThePEG/Utilities/Debug.h"
#include "Herwig/Shower/QTilde/QTildeShowerHandler.h"
#include "Herwig/Shower/QTilde/Base/PartnerFinder.h"
#include "Herwig/Shower/QTilde/Kinematics/KinematicsReconstructor.h"
#include "Herwig/Shower/QTilde/Base/ShowerVertex.h"
#include <cassert>

using namespace Herwig;

void IS_QTildeShowerKinematics1to2::
updateChildren( const tShowerParticlePtr theParent, 
		const ShowerParticleVector & children,
		ShowerPartnerType) const {
  const ShowerParticle::Parameters & parent = theParent->showerParameters();
  ShowerParticle::Parameters & child0 = children[0]->showerParameters();
  ShowerParticle::Parameters & child1 = children[1]->showerParameters();
  double cphi = cos(phi());
  double sphi = sin(phi());

  child1.alpha = (1.-z()) * parent.alpha;

  child1.ptx = (1.-z()) * parent.ptx - cphi * pT();
  child1.pty = (1.-z()) * parent.pty - sphi * pT();
  child1.pt  = sqrt( sqr(child1.ptx) + sqr(child1.pty) );
  // space-like child
  child0.alpha = parent.alpha - child1.alpha;
  child0.beta  = parent.beta  - child1.beta;
  child0.ptx   = parent.ptx   - child1.ptx;
  child0.pty   = parent.pty   - child1.pty;
}


void IS_QTildeShowerKinematics1to2::
updateParent(const tShowerParticlePtr parent, 
	     const ShowerParticleVector & children,
	     unsigned int ,
	     ShowerPartnerType partnerType) const {
  // calculate the scales
  splittingFn()->evaluateInitialStateScales(partnerType,scale(),z(),parent,
					    children[0],children[1]);
  // set proper colour connections
  splittingFn()->colourConnection(parent,children[0],children[1],
				  partnerType,true);
  // set proper parent/child relationships
  parent->addChild(children[0]);
  parent->addChild(children[1]);
  parent->x(children[0]->x()/z());
  // sort out the helicity stuff 
  // construct the spin info for parent and timelike child
  // temporary assignment of shower parameters to calculate correlations
  parent->showerParameters().alpha = parent->x();
  children[1]->showerParameters().alpha = (1.-z()) * parent->x();
  children[1]->showerParameters().ptx   = - cos(phi()) * pT();
  children[1]->showerParameters().pty   = - sin(phi()) * pT();
  children[1]->showerParameters().pt    = pT();
  parent     ->showerBasis(children[0]->showerBasis(),true);
  children[1]->showerBasis(children[0]->showerBasis(),true);
  parent     ->setShowerMomentum(false);
  children[1]->setShowerMomentum(true);
  if(! dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->correlations()) return;
  SpinPtr pspin(children[0]->spinInfo());
  if(!pspin ||  !dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->spinCorrelations() ) return;
  // compute the matrix element for spin correlations
  IdList ids;
  ids.push_back(parent->dataPtr());
  ids.push_back(children[0]->dataPtr());
  ids.push_back(children[1]->dataPtr());
  Energy2 t = (1.-z())*sqr(scale())/z();
  // create the vertex
  SVertexPtr vertex(new_ptr(ShowerVertex()));
  // set the matrix element
  vertex->ME(splittingFn()->matrixElement(z(),t,ids,phi(),false));
  // set the incoming particle for the vertex 
  // (in reality the first child as going backwards)
  pspin->decayVertex(vertex);
  // construct the spin infos
  parent     ->constructSpinInfo(false);
  children[1]->constructSpinInfo(true);
  // connect the spinInfo objects to the vertex
  parent     ->spinInfo()->productionVertex(vertex);
  children[1]->spinInfo()->productionVertex(vertex);
}

void IS_QTildeShowerKinematics1to2::
reconstructParent(const tShowerParticlePtr parent, 
		  const ParticleVector & children ) const {
  PPtr c1 = children[0];
  ShowerParticlePtr c2 = dynamic_ptr_cast<ShowerParticlePtr>(children[1]);
  ShowerParticle::Parameters & c2param = c2->showerParameters();
  // get shower variables from 1st child in order to keep notation
  // parent->(c1, c2) clean even though the splitting was initiated
  // from c1.  The name updateParent is still referring to the
  // timelike branching though.
  // on-shell child
  
  auto m=  ShowerHandler::currentHandler()->retConstituentMasses()?
           c2->data().constituentMass():
           c2->data().mass();
  
  c2param.beta = 0.5*( sqr(m) + sqr(c2param.pt) )
    / ( c2param.alpha * parent->showerBasis()->p_dot_n() );
  Lorentz5Momentum pnew = parent->showerBasis()->
    sudakov2Momentum(c2param.alpha, c2param.beta, 
		     c2param.ptx  , c2param.pty);
  pnew.setMass(m);
  pnew.rescaleEnergy();
  c2->set5Momentum( pnew );
  // spacelike child
  Lorentz5Momentum pc1(parent->momentum() - c2->momentum());
  pc1.rescaleMass();
  c1->set5Momentum(pc1);
}

void IS_QTildeShowerKinematics1to2::
updateLast( const tShowerParticlePtr theLast,Energy px,Energy py) const {
  if(theLast->isFinalState()) return;
  ShowerParticle::Parameters & last = theLast->showerParameters();
  Lorentz5Momentum pVector = theLast->showerBasis()->pVector();
  Energy2 pt2 = sqr(px) + sqr(py);
  last.alpha = theLast->x();
  last.beta  = 0.5 * pt2 / last.alpha / theLast->showerBasis()->p_dot_n();
  last.ptx   = ZERO;
  last.pty   = ZERO;
  last.pt    = ZERO;
  // momentum
  Lorentz5Momentum ntemp = Lorentz5Momentum(ZERO,-pVector.vect());
  double beta = 0.5 * pt2 / last.alpha / ( pVector * ntemp);
  Lorentz5Momentum plast = 
    Lorentz5Momentum( (pVector.z()>ZERO ? px : -px), py, ZERO, ZERO)
    + theLast->x() * pVector + beta * ntemp;
  plast.rescaleMass();
  theLast->set5Momentum(plast);
}
