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
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/LorentzSpinorBar.h"
#include "Herwig++/Shower/ShowerHandler.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/ShowerModel.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

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
  // sort out the helicity stuff 
  if(! ShowerHandler::currentHandler()->evolver()->correlations()) return;
  SpinPtr pspin(parent->spinInfo());
  if(!pspin) return;
  // get the vertex
  VertexPtr vertex(const_ptr_cast<VertexPtr>(pspin->decayVertex()));
  if(!vertex) return;
  // construct the spin info for the children
  ShowerParticleVector::const_iterator pit;
  for(pit=children.begin();pit!=children.end();++pit) {
    Energy mass = (*pit)->data().mass(); 
    // calculate the momentum of the children assuming on-shell
    Energy2 pt2 = sqr((**pit).showerParameters().pt);
    double alpha = (**pit).showerParameters().alpha;
    double beta = 0.5*(sqr(mass) + pt2 - sqr(alpha)*pVector().m2())/(alpha*p_dot_n());
    Lorentz5Momentum porig=sudakov2Momentum(alpha,beta,
					    (**pit).showerParameters().ptx,
					    (**pit).showerParameters().pty);
    porig.setMass(mass);
    // now construct the required spininfo and calculate the basis states
    PDT::Spin spin((*pit)->dataPtr()->iSpin());
    if(spin==PDT::Spin0) {
      assert(false);
    }
    // calculate the basis states and construct the SpinInfo for a spin-1/2 particle
    else if(spin==PDT::Spin1Half) {
      // outgoing particle
      if((*pit)->id()>0) {
      	vector<LorentzSpinorBar<SqrtEnergy> > stemp;
	SpinorBarWaveFunction::calculateWaveFunctions(stemp,*pit,outgoing);
	SpinorBarWaveFunction::constructSpinInfo(stemp,*pit,outgoing,true);
      }
      // outgoing antiparticle
      else {
      	vector<LorentzSpinor<SqrtEnergy> > stemp;
	SpinorWaveFunction::calculateWaveFunctions(stemp,*pit,outgoing);
	SpinorWaveFunction::constructSpinInfo(stemp,*pit,outgoing,true);
      }
    }
    // calculate the basis states and construct the SpinInfo for a spin-1 particle
    else if(spin==PDT::Spin1) {
      bool massless((*pit)->id()==ParticleID::g||(*pit)->id()==ParticleID::gamma);
      vector<Helicity::LorentzPolarizationVector> vtemp;
      VectorWaveFunction::calculateWaveFunctions(vtemp,*pit,outgoing,massless);
      VectorWaveFunction::constructSpinInfo(vtemp,*pit,outgoing,true,massless);
    }
    else {
      throw Exception() << "Spins higher than 1 are not yet implemented in " 
			<< "FS_QtildaShowerKinematics1to2::constructVertex() "
			<< Exception::runerror;
    }
    // connect the spinInfo object to the vertex
    (*pit)->spinInfo()->productionVertex(vertex);
    (*pit)->set5Momentum(porig);
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
  parent->set5Momentum( c1->momentum() + c2->momentum() );
}

void FS_QTildeShowerKinematics1to2::reconstructLast(const tShowerParticlePtr theLast,
						    Energy mass) const {
  // set beta component and consequently all missing data from that,
  // using the nominal (i.e. PDT) mass.
  Energy theMass = mass > ZERO  ?  mass : theLast->data().constituentMass();
  ShowerParticle::Parameters & last = theLast->showerParameters();
  last.beta = ( sqr(theMass) + sqr(last.pt) - sqr(last.alpha) * pVector().m2() )
    / ( 2. * last.alpha * p_dot_n() );
  // set that new momentum
  theLast->set5Momentum(sudakov2Momentum( last.alpha, last.beta, 
					  last.ptx, last.pty) );
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
  Energy2 q2 = 
    sqr(children[0]->virtualMass())/z() + 
    sqr(children[1]->virtualMass())/(1.-z()) +
    pt2/z()/(1.-z());
  if(pt2<ZERO) {
    parent->virtualMass(ZERO);
  }
  else {
    parent->virtualMass(sqrt(q2));
    pT(sqrt(pt2));
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
