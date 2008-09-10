// -*- C++ -*-
//
// IS_QtildaShowerKinematics1to2.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_QtildaShowerKinematics1to2 class.
//

#include "IS_QtildaShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

void IS_QtildaShowerKinematics1to2::
updateChildren( const tShowerParticlePtr theParent, 
		const ShowerParticleVector theChildren ) const {
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


void IS_QtildaShowerKinematics1to2::
updateParent(const tShowerParticlePtr theParent, 
	     const ShowerParticleVector theChildren ) const {
  // no z for angular ordering in backward branchings
  theParent->setEvolutionScale(scale());
  theChildren[1]->setEvolutionScale((1.-z())*scale());
  // set proper colour connections
  splittingFn()->colourConnection(theParent,theChildren[0],theChildren[1],true);
  // set proper parent/child relationships
  theParent->addChild(theChildren[0]);
  theParent->addChild(theChildren[1]);
  theParent->x(theChildren[0]->x()/z());
  // create the storage for the shower variables
  theParent->showerVariables().resize(3);
  theParent->showerParameters().resize(2);
  if(theChildren[0]->showerVariables().empty()) {
    theChildren[0]->showerVariables().resize(3);
    theChildren[0]->showerParameters().resize(2);
  }
}

void IS_QtildaShowerKinematics1to2::
reconstructParent(const tShowerParticlePtr theParent, 
		  const ParticleVector theChildren ) const {
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

void IS_QtildaShowerKinematics1to2::
updateLast( const tShowerParticlePtr theLast,Energy px,Energy py) const {
  if(theLast->isFinalState()) return;
  Energy2 pt2=sqr(px)+sqr(py);
  theLast->showerParameters()[0]=theLast->x();
  theLast->showerParameters()[1]=0.5*pt2/
    theLast->showerParameters()[0]/p_dot_n();
  theLast->showerVariables().resize(3);
  theLast->showerParameters().resize(2);
  for(unsigned int ix=0;ix<3;++ix) theLast->showerVariables()[ix]=0.*MeV;
  // momentum
  Lorentz5Momentum ntemp=Lorentz5Momentum(0.*GeV,-pVector().vect());
  double beta = 0.5*pt2/
    theLast->showerParameters()[0]/(pVector()*ntemp);
  Lorentz5Momentum plast = 
    Lorentz5Momentum(pVector().z()>0.*GeV ? px : -px ,py,0.*GeV,0.*GeV)
    +theLast->x()*pVector()+beta*ntemp;
  plast.rescaleMass();
  theLast->set5Momentum(plast);
}
 
void IS_QtildaShowerKinematics1to2::initialize(ShowerParticle & particle, PPtr parent) {
  // For the time being we are considering only 1->2 branching
  Lorentz5Momentum p, n, pthis, ppartner, pcm;
  assert(particle.perturbative()!=2);
  if(particle.perturbative()==1) {
    // find the partner and its momentum
    ShowerParticlePtr partner=particle.partner();
    Lorentz5Momentum ppartner(partner->momentum());
    if(partner->getThePEGBase()) ppartner=partner->getThePEGBase()->momentum();
    if(partner->isFinalState()) {
      Lorentz5Momentum pa = -particle.momentum()+partner->momentum();
      Lorentz5Momentum pb =  particle.momentum();
      Axis axis(pa.vect().unit());
      LorentzRotation rot;
      double sinth(sqrt(1.-sqr(axis.z())));
      rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      rot.rotateX(Constants::pi);
      rot.boostZ( pa.e()/pa.vect().mag());
      pb*=rot;
      Boost trans = -1./pb.e()*pb.vect();
      trans.setZ(0.);
      rot.boost(trans);
      pcm = rot*parent->momentum();
      rot.invert();
      n = rot*Lorentz5Momentum(0.*GeV,-pcm.vect());
      p = rot*Lorentz5Momentum(0.*GeV, pcm.vect());
    }
    else {
      pcm = parent->momentum();
      p = Lorentz5Momentum(0.0*MeV, pcm.vect());
      n = Lorentz5Momentum(0.0*MeV, -pcm.vect());
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
