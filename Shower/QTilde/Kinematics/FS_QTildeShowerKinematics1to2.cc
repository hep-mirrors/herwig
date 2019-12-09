// -*- C++ -*-
//
// FS_QTildeShowerKinematics1to2.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_QTildeShowerKinematics1to2 class.
//

#include "FS_QTildeShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Shower/QTilde/SplittingFunctions/SplittingFunction.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"
#include "ThePEG/Utilities/Debug.h"
#include "Herwig/Shower/QTilde/QTildeShowerHandler.h"
#include "Herwig/Shower/QTilde/Base/PartnerFinder.h"
#include "Herwig/Shower/QTilde/Kinematics/KinematicsReconstructor.h"
#include "Herwig/Shower/QTilde/Kinematics/KinematicHelpers.h"
#include "Herwig/Shower/QTilde/Base/ShowerVertex.h"

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
	       ShowerPartnerType partnerType) const {
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
    (**pit).showerBasis(parent->showerBasis(),true);
    (**pit).setShowerMomentum(true);
  }
  // sort out the helicity stuff 
  if(! dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->correlations()) return;
  SpinPtr pspin(parent->spinInfo());
  if(!pspin ||  !dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->spinCorrelations() ) return;
  Energy2 t = sqr(scale())*z()*(1.-z());
  IdList ids;
  ids.push_back(parent->dataPtr());
  ids.push_back(children[0]->dataPtr());
  ids.push_back(children[1]->dataPtr());
  // create the vertex
  SVertexPtr vertex(new_ptr(ShowerVertex()));
  // set the matrix element
  vertex->ME(splittingFn()->matrixElement(z(),t,ids,phi(),true));
  RhoDMatrix mapping;
  SpinPtr inspin;
  bool needMapping = parent->getMapping(inspin,mapping);
  if(needMapping) vertex->incomingBasisTransform(mapping);
  // set the incoming particle for the vertex
  parent->spinInfo()->decayVertex(vertex);
  for(ShowerParticleVector::const_iterator pit=children.begin();
      pit!=children.end();++pit) {
    // construct the spin info for the children
    (**pit).constructSpinInfo(true);
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
  Energy theMass =ZERO;
  if(!(mass > ZERO) && ShowerHandler::currentHandler()->retConstituentMasses())
    theMass = last->data().constituentMass();
  else
    theMass = mass > ZERO ? mass : last->data().mass();
  
  Lorentz5Momentum pVector = last->showerBasis()->pVector();
  ShowerParticle::Parameters & lastParam = last->showerParameters();
  Energy2 denom = 2. * lastParam.alpha * last->showerBasis()->p_dot_n();
  if(abs(denom)/(sqr(pVector.e())+pVector.rho2())<1e-10) {
    throw KinematicsReconstructionVeto();
  }
  lastParam.beta = ( sqr(theMass) + sqr(lastParam.pt)
		     - sqr(lastParam.alpha) * pVector.m2() )
    / denom;
  // set that new momentum
  Lorentz5Momentum newMomentum = last->showerBasis()->
    sudakov2Momentum( lastParam.alpha, lastParam.beta,
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

void FS_QTildeShowerKinematics1to2::updateParent(const tShowerParticlePtr parent, 
						 const ShowerParticleVector & children,
						 unsigned int pTscheme,
						 ShowerPartnerType) const {
  IdList ids(3);
  ids[0] = parent->dataPtr();
  ids[1] = children[0]->dataPtr();
  ids[2] = children[1]->dataPtr();
  const vector<Energy> & virtualMasses = SudakovFormFactor()->virtualMasses(ids);
  if(children[0]->children().empty()) children[0]->virtualMass(virtualMasses[1]);
  if(children[1]->children().empty()) children[1]->virtualMass(virtualMasses[2]);
  // compute the new pT of the branching
  Energy2 m02 = (ids[0]->id()!=ParticleID::g && ids[0]->id()!=ParticleID::gamma) ?
  	sqr(virtualMasses[0]) : Energy2();
  Energy2 pt2;
  if(pTscheme==0) {
    pt2 = QTildeKinematics::pT2_FSR(sqr(scale()), z(), m02,
				    sqr(virtualMasses[1])          ,sqr(virtualMasses[2]),
				    sqr(virtualMasses[1])          ,sqr(virtualMasses[2]));
  }
  else if(pTscheme==1) {
    pt2 = QTildeKinematics::pT2_FSR(sqr(scale()), z(), m02,
				    sqr(virtualMasses[1])          ,sqr(virtualMasses[2]),
				    sqr(children[0]->virtualMass()), sqr(children[1]->virtualMass()));
  }
  else if(pTscheme==2) {
    pt2 = QTildeKinematics::pT2_FSR(sqr(scale()), z(), m02,
				    sqr(children[0]->virtualMass()), sqr(children[1]->virtualMass()),
				    sqr(children[0]->virtualMass()), sqr(children[1]->virtualMass()));
  }
  else
    assert(false);
  
  if(pt2>ZERO) {
    pT(sqrt(pt2));
  }
  else {
    pt2=ZERO;
    pT(ZERO);
  }
  Energy2 q2 = QTildeKinematics::q2_FSR(
      pt2, z(), sqr(children[0]->virtualMass()), sqr(children[1]->virtualMass())
    );
  parent->virtualMass(sqrt(q2));
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
