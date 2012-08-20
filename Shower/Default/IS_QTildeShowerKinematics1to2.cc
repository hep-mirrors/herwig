// -*- C++ -*-
//
// IS_QTildeShowerKinematics1to2.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_QTildeShowerKinematics1to2 class.
//

#include "IS_QTildeShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/Utilities/Debug.h"
#include <cassert>

using namespace Herwig;

void IS_QTildeShowerKinematics1to2::
updateChildren( const tShowerParticlePtr theParent, 
		const ShowerParticleVector & theChildren,
		ShowerPartnerType::Type) const {
  theChildren[1]->showerVariables().resize(3);
  theChildren[1]->showerParameters().resize(2);
  theChildren[1]->showerParameters()[0]=
    (1.-z())*theParent->showerParameters()[0];
  double cphi = cos(phi()),sphi = sin(phi());
  theChildren[1]->showerVariables()[0] = 
    (1.-z())*theParent->showerVariables()[0] - cphi*pT();
  theChildren[1]->showerVariables()[1]=
    (1.-z())*theParent->showerVariables()[1] - sphi*pT();
  theChildren[1]->showerVariables()[2]=
    sqrt(sqr(theChildren[1]->showerVariables()[0])+
	 sqr(theChildren[1]->showerVariables()[1]));
  // space-like child
  theChildren[0]->showerParameters()[0] =
    theParent->showerParameters()[0] - theChildren[1]->showerParameters()[0];
  theChildren[0]->showerParameters()[1] = 
    theParent->showerParameters()[1] - theChildren[1]->showerParameters()[1];
  theChildren[0]->showerVariables()[0]=
    theParent->showerVariables()[0] - theChildren[1]->showerVariables()[0];
  theChildren[0]->showerVariables()[1]=
    theParent->showerVariables()[1] - theChildren[1]->showerVariables()[1];
}


void IS_QTildeShowerKinematics1to2::
updateParent(const tShowerParticlePtr parent, 
	     const ShowerParticleVector & children,
	     ShowerPartnerType::Type partnerType) const {
  // calculate the scales
  splittingFn()->evaluateInitialStateScales(partnerType,scale(),z(),parent,
					    children[0],children[1]);
  // debugging printout if needed
  if(Debug::level >= 10 ) printScales(parent,children[0],children[1]);
  // set proper colour connections
  splittingFn()->colourConnection(parent,children[0],children[1],
				  partnerType,true);
  // set proper parent/child relationships
  parent->addChild(children[0]);
  parent->addChild(children[1]);
  parent->x(children[0]->x()/z());
  // create the storage for the shower variables
  parent->showerVariables().resize(3);
  parent->showerParameters().resize(2);
  if(children[0]->showerVariables().empty()) {
    children[0]->showerVariables().resize(3);
    children[0]->showerParameters().resize(2);
  }
}

void IS_QTildeShowerKinematics1to2::
reconstructParent(const tShowerParticlePtr theParent, 
		  const ParticleVector & theChildren ) const {
  ShowerParticlePtr c1 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[0]);
  ShowerParticlePtr c2 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[1]);
  // get shower variables from 1st child in order to keep notation
  // parent->(c1, c2) clean even though the splitting was initiated
  // from c1.  The name updateParent is still referring to the
  // timelike branching though.
  // on-shell child
  c2->showerParameters()[1]=
    0.5*(sqr(c2->data().constituentMass())+sqr(c2->showerVariables()[2]))
    /(c2->showerParameters()[0]*p_dot_n());
  c2->set5Momentum(sudakov2Momentum(c2->showerParameters()[0],
				    c2->showerParameters()[1], 
				    c2->showerVariables()[0],
				    c2->showerVariables()[1],0));
  // spacelike child
  Lorentz5Momentum pc1(theParent->momentum() - c2->momentum());
  pc1.rescaleMass();
  c1->set5Momentum(pc1);
}

void IS_QTildeShowerKinematics1to2::
updateLast( const tShowerParticlePtr theLast,Energy px,Energy py) const {
  if(theLast->isFinalState()) return;
  Energy2 pt2=sqr(px)+sqr(py);
  theLast->showerParameters()[0]=theLast->x();
  theLast->showerParameters()[1]=0.5*pt2/
    theLast->showerParameters()[0]/p_dot_n();
  theLast->showerVariables().resize(3);
  theLast->showerParameters().resize(2);
  for(unsigned int ix=0;ix<3;++ix) theLast->showerVariables()[ix]=ZERO;
  // momentum
  Lorentz5Momentum ntemp=Lorentz5Momentum(ZERO,-pVector().vect());
  double beta = 0.5*pt2/
    theLast->showerParameters()[0]/(pVector()*ntemp);
  Lorentz5Momentum plast = 
    Lorentz5Momentum(pVector().z()>ZERO ? px : -px ,py,ZERO,ZERO)
    +theLast->x()*pVector()+beta*ntemp;
  plast.rescaleMass();
  theLast->set5Momentum(plast);
}
 
void IS_QTildeShowerKinematics1to2::initialize(ShowerParticle & particle, PPtr parent) {
  // For the time being we are considering only 1->2 branching
  Lorentz5Momentum p, n, pthis, pcm;
  assert(particle.perturbative()!=2);
  if(particle.perturbative()==1) {
    // find the partner and its momentum
    ShowerParticlePtr partner=particle.partner();
    assert(partner);
    if(partner->isFinalState()) {
      Lorentz5Momentum pa = -particle.momentum()+partner->momentum();
      Lorentz5Momentum pb =  particle.momentum();
      Energy scale=parent->momentum().t();
      Lorentz5Momentum pbasis(ZERO,parent->momentum().vect().unit()*scale);
      Axis axis(pa.vect().unit());
      LorentzRotation rot;
      double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
      if(axis.perp2()>1e-20) {
	rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
	rot.rotateX(Constants::pi);
      }
      if(abs(1.-pa.e()/pa.vect().mag())>1e-6) rot.boostZ( pa.e()/pa.vect().mag());
      pb *= rot;


      if(pb.perp2()/GeV2>1e-20) {
	Boost trans = -1./pb.e()*pb.vect();
	trans.setZ(0.);
	rot.boost(trans);
      }
      pbasis *=rot;
      rot.invert();
      n = rot*Lorentz5Momentum(ZERO,-pbasis.vect());
      p = rot*Lorentz5Momentum(ZERO, pbasis.vect());
    }
    else {
      pcm = parent->momentum();
      p = Lorentz5Momentum(ZERO, pcm.vect());
      n = Lorentz5Momentum(ZERO, -pcm.vect());
    }
  } 
  else {
    p = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
      ->showerKinematics()->getBasis()[0];
    n = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
      ->showerKinematics()->getBasis()[1];
  }
  setBasis(p,n);
}
