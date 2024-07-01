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
#include "Herwig/Shower/QTilde/Kinematics/KinematicHelpers.h"
#include "Herwig/Shower/QTilde/Base/ShowerVertex.h"
#include "Herwig/Shower/QTilde/SplittingFunctions/Sudakov1to2FormFactor.h"
#include <cassert>

using namespace Herwig;

IS_QTildeShowerKinematics1to2::IS_QTildeShowerKinematics1to2(Energy scale, double z, double phi, Energy pt, tSudakovPtr sud) 
  : ShowerKinematics(scale,z,phi,pt,sud), sudakov1to2_(dynamic_ptr_cast<tSudakov1to2Ptr>(sud)) {}

void IS_QTildeShowerKinematics1to2::
updateChildren( const tShowerParticlePtr theParent, 
		const ShowerParticleVector & children,
		unsigned int pTscheme,
		ShowerPartnerType) const {
  const ShowerParticle::Parameters & parent = theParent->showerParameters();
  ShowerParticle::Parameters & child0 = children[0]->showerParameters();
  ShowerParticle::Parameters & child1 = children[1]->showerParameters();
  double cphi = cos(phi());
  double sphi = sin(phi());
  ///-- Perform the reconstruction depenting on the recoil scheme, following 2107.04051 (Bewick, Ferrario Ravasio, Richardson and Seymour)
  Energy2 m02 = ZERO;
  Energy2 m12 = ZERO;
  // Time like children, which can be massive
  IdList ids(1);
  ids[0] =  children[1]->dataPtr();
  Energy2 m22 = ZERO;
  if (ids[0]->id()!=ParticleID::g && ids[0]->id()!=ParticleID::gamma) {
    const vector<Energy> & onshellMasses = sudakov1to2_->virtualMasses(ids);
    m22 =  sqr(onshellMasses[0]);
  }

  
  if(theParent->parents().empty()) theParent->virtualMass(ZERO);
  if(children[1]->children().empty()) children[1]->virtualMass(sqrt(m22));


  // each particle contains a pointer to its mass, that we can use to calculate the virtuality.
  // onshell incoming partons are massless, while offshell ones have negative virtuality q2
  // q2 = - virtualMass**2  for ISR [THIS IS NOT TRUE FOR RESONANCES]
  Energy2 pt2;
  if(pTscheme==0) {
    pt2 = QTildeKinematics::pT2_ISR_new(sqr(scale()), z(), m02,
					m12          ,m22,
					m02          ,m22);
  }
  else if(pTscheme==1) {
    pt2 = QTildeKinematics::pT2_ISR_new(sqr(scale()), z(), m02,
					m12          ,m22,
					-sqr(theParent->virtualMass()), sqr(children[1]->virtualMass()));
  }
  else if(pTscheme==2) {
    pt2 = QTildeKinematics::pT2_ISR_new(sqr(scale()), z(), -sqr(theParent->virtualMass()),
					m12          ,sqr(children[1]->virtualMass()),
				        -sqr(theParent->virtualMass()) ,sqr(children[1]->virtualMass()));
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
  Energy2 q2 = QTildeKinematics::q2_ISR_new(pt2, z(), -sqr(theParent->virtualMass()) ,sqr(children[1]->virtualMass()));
  // the virtuality is negative
  if(q2/sqr(1_GeV) > 0.){
      std::cout<<" virtuality of spacelike parton >0 !! ERROR!! Q2="<<q2/sqr(1_GeV)<<std::endl;
      abort();
    }
  children[0]->virtualMass(sqrt(-q2));
  // -- now as before 2107.04051

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
	     unsigned int,
	     ShowerPartnerType partnerType) const {
  // calculate the scales
  sudakov1to2_->evaluateInitialStateScales(partnerType,scale(),z(),parent,
					    children[0],children[1]);
  // set proper colour connections
  sudakov1to2_->colourConnection(parent,children[0],children[1],
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
  vertex->ME(sudakov1to2_->matrixElement(z(),t,ids,phi(),false));
  // set the incoming particle for the vertex 
  // (in reality the first child as going backwards)
  pspin->decayVertex(vertex);
  // check if rotating basis
  RhoDMatrix mapping;
  SpinPtr output;
  bool needMapping = children[0]->getMapping(output,mapping);
  if(needMapping)
    vertex->outgoingBasisTransform(0,mapping);
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
  // fix the spin info AHHHH
  if(c2->spinInfo()&&abs(c2->id())==ParticleID::tauminus) {
    // rotate old momentum along z
    Axis axis(c2->momentum().vect().unit());
    double sinth(sqrt(1.-sqr(axis.z())));
    LorentzRotation r1;
    if(axis.perp2()>1e-10) {
      r1.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    }
    Lorentz5Momentum p1 = r1*c2->momentum();
    // rotate new momentum along z
    axis = pnew.vect().unit();
    sinth = sqrt(1.-sqr(axis.z()));
    LorentzRotation r2;
    if(axis.perp2()>1e-10) {
      r2.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    }
    Lorentz5Momentum p2 = r2*pnew;
    // compute boost along z old -> new
    // full result
    double a = p2.t()/p2.z();
    double b1  = (p1.t()-a*p1.z())/(p1.z()-a*p1.t());
    // small mass expansion
    Energy2 o2 = sqr(p1.z()), n2=sqr(p2.z());
    double b2  = (o2-n2)/(o2+n2)*(1. - sqr(p1.mass())/(n2+o2) +pow<4,1>(p1.mass())*(0.125/(n2*o2)+1./sqr(n2+o2)));
    // compute the total boost
    LorentzRotation r3; r3.boostZ(b1);
    LorentzRotation r4; r4.boostZ(b2);
    Energy t1 = p2.t(), t2 = (r3*p1).t(), t3 = (r4*p1).t();
    LorentzRotation r = abs(t1-t2)<abs(t1-t3) ? r3*r1 : r4*r1;
    r *= r2.inverse();
    c2->transform(r);
  }
  else
    c2->set5Momentum( pnew );
  // spacelike child
  Lorentz5Momentum pc1(parent->momentum() - c2->momentum());
  pc1.rescaleMass();
  c1->set5Momentum(pc1);
}

void IS_QTildeShowerKinematics1to2::
updateLast( const tShowerParticlePtr theLast,Energy px,Energy py) const {
  if(theLast->isFinalState()) return; 
  Lorentz5Momentum pVector = theLast->showerBasis()->pVector();
  ShowerParticle::Parameters & last = theLast->showerParameters();
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


void IS_QTildeShowerKinematics1to2::
resetChildren(const tShowerParticlePtr parent, const ShowerParticleVector &) const {

  if(parent->children().size() == 2){

    if(parent->children()[1]->children().size()==2){
       tShowerParticlePtr timelikeChild = dynamic_ptr_cast<ShowerParticlePtr>(parent->children()[1]);
       ShowerParticleVector timelikeGrandChildren;
       for(unsigned int iy=0;iy<parent->children()[1]->children().size();++iy)
	 timelikeGrandChildren.push_back(dynamic_ptr_cast<ShowerParticlePtr>
       					 (parent->children()[1]->children()[iy]));
       //timelike reset children
       timelikeChild->showerKinematics()->resetChildren(timelikeChild, timelikeGrandChildren);
    }
    
    if(parent->children()[0]->children().size()==2){
      const tShowerParticlePtr spacelikeChild = dynamic_ptr_cast<ShowerParticlePtr>(parent->children()[0]);
      //Check if its timelike child is a valid shower particle istance. If not, we end here the reshuffling.
      tShowerParticlePtr timelikeChild = dynamic_ptr_cast<ShowerParticlePtr>(spacelikeChild->children()[1]);
      if(timelikeChild){
      ShowerParticleVector dummyChildren;
      spacelikeChild->showerKinematics()->resetChildren(spacelikeChild, dummyChildren);
      }
    }
          
  }
  
}
	
