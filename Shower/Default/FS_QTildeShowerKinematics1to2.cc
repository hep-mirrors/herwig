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

using namespace Herwig;

void FS_QTildeShowerKinematics1to2::
updateChildren(const tShowerParticlePtr parent, 
	       const ShowerParticleVector & children,
	       ShowerPartnerType::Type partnerType, 
	       bool angularOrder) const {
  assert(children.size()==2);
  // resize the parameter vectors
  if(parent->showerVariables().empty()) {
    parent->showerVariables().resize(3);
    parent->showerParameters().resize(2);
    parent->showerParameters()[0]=1.;
  }
  children[0]->showerVariables() .resize(3);
  children[0]->showerParameters().resize(2);
  children[1]->showerVariables() .resize(3);
  children[1]->showerParameters().resize(2);
  // identify emittor and spectator
  tShowerParticlePtr emittor = children[0];
  tShowerParticlePtr emitted = children[1];
  double zChild[2] = {z(),1.-z()};
  bool bosonSplitting(false);
  // special for g -> gg, particle highest z is emittor
  if(emittor->id()==emitted->id()&&emittor->id()==parent->id()) {
    if(zChild[1]>zChild[0]) {
      swap(zChild[0],zChild[1]);
      swap(emitted,emittor);
    }
  }
  // otherwise if particle ID same
  else if(emitted->id()==parent->id()) {
    swap(zChild[0],zChild[1]);
    swap(emitted,emittor);
  }
  // no real emittor/eemitted
  else if(emittor->id()!=parent->id()) {
    bosonSplitting = true;
  }
  // scales for angular ordering
  Energy AOScale[2];
  if(angularOrder) {
    for(unsigned int ix=0;ix<2;++ix) AOScale [ix] = zChild[ix]*scale(); 
  }
  else {
    for(unsigned int ix=0;ix<2;++ix) AOScale [ix] =            scale();
  }
  // if not these cases doesn't matter
  // now the various scales
  if(partnerType==ShowerPartnerType::QED) {
    // normal case
    if(!bosonSplitting) {
      for(map<ShowerPartnerType::Type,pair<Energy,Energy> >::const_iterator 
	    it = parent->evolutionScales().begin();
	  it!=parent->evolutionScales().end();++it) {
	emittor->evolutionScale(it->first,make_pair(min(AOScale[0],it->second.first ),
						    min(scale()   ,it->second.second)));
	if(it->first==ShowerPartnerType::QED) {
	  emitted->evolutionScale(it->first,make_pair(min(AOScale[1],it->second.first ),
						      min(scale()   ,it->second.second)));
	}
      }
    }
    // gamma -> f fbar
    else {
      // set QED scales as normal
      for(map<ShowerPartnerType::Type,pair<Energy,Energy> >::const_iterator 
	    it = parent->evolutionScales().begin();
	  it!=parent->evolutionScales().end();++it) {
	assert(it->first==ShowerPartnerType::QED);
	emittor->evolutionScale(it->first,make_pair(AOScale[0],scale()));
	emitted->evolutionScale(it->first,make_pair(AOScale[1],scale()));
      }
      // and any QCD scales needed
      PDT::Colour emittorColour = emittor->dataPtr()->iColour();
      if(emittorColour==PDT::Colour3||emittorColour==PDT::Colour8||emittorColour==PDT::Colour6) {
	emittor->evolutionScale(ShowerPartnerType::QCDColourLine,
				make_pair(AOScale[0],scale()));
      }
      if(emittorColour==PDT::Colour3bar||emittorColour==PDT::Colour8||emittorColour==PDT::Colour6bar) {
	emittor->evolutionScale(ShowerPartnerType::QCDAntiColourLine,
				make_pair(AOScale[0],scale()));
      }
      PDT::Colour emittedColour = emitted->dataPtr()->iColour();
      if(emittedColour==PDT::Colour3||emittedColour==PDT::Colour8||emittedColour==PDT::Colour6) {
	emitted->evolutionScale(ShowerPartnerType::QCDColourLine,
				make_pair(AOScale[1],scale()));
      }
      if(emittedColour==PDT::Colour3bar||emittedColour==PDT::Colour8||emittedColour==PDT::Colour6bar) {
	emitted->evolutionScale(ShowerPartnerType::QCDAntiColourLine,
				make_pair(AOScale[1],scale()));
      }
    }
  }
  else {
    // normal case
    if(!bosonSplitting) {
      // scales for the emittor
      for(map<ShowerPartnerType::Type,pair<Energy,Energy> >::const_iterator 
	    it = parent->evolutionScales().begin();
	  it!=parent->evolutionScales().end();++it) {
	emittor->evolutionScale(it->first,make_pair(min(AOScale[0],it->second.first ),
						    min(scale()   ,it->second.second)));
      }
    }
    else {
      // QCD scales needed
      PDT::Colour emittorColour = emittor->dataPtr()->iColour();
      if(emittorColour==PDT::Colour3||emittorColour==PDT::Colour8||emittorColour==PDT::Colour6) {
	emittor->evolutionScale(ShowerPartnerType::QCDColourLine,
				make_pair(AOScale[0],scale()));
      }
      if(emittorColour==PDT::Colour3bar||emittorColour==PDT::Colour8||emittorColour==PDT::Colour6bar) {
	emittor->evolutionScale(ShowerPartnerType::QCDAntiColourLine,
				make_pair(AOScale[0],scale()));
      }
      // QED scales
      if(emittor->dataPtr()->charged()) {
	emittor->evolutionScale(ShowerPartnerType::QED,make_pair(AOScale[0],scale()));
      }
    }
    PDT::Colour emittedColour = emitted->dataPtr()->iColour();
    if(emittedColour==PDT::Colour3||emittedColour==PDT::Colour8||emittedColour==PDT::Colour6) {
      emitted->evolutionScale(ShowerPartnerType::QCDColourLine,
			      make_pair(AOScale[1],scale()));
    }
    if(emittedColour==PDT::Colour3bar||emittedColour==PDT::Colour8||emittedColour==PDT::Colour6bar) {
      emitted->evolutionScale(ShowerPartnerType::QCDAntiColourLine,
			      make_pair(AOScale[1],scale()));
    }
    if(emitted->dataPtr()->charged()) {
      emitted->evolutionScale(ShowerPartnerType::QED,make_pair(AOScale[1],scale()));
    }
  }
  // debugging printout ifd needed
  if(Debug::level >= 10 ) printScales(parent,children[0],children[1]);
  // determine alphas of children according to interpretation of z
  children[0]->showerParameters()[0]=     z() *parent->showerParameters()[0];
  children[1]->showerParameters()[0]= (1.-z())*parent->showerParameters()[0];
  // set the values
  children[0]->showerVariables()[0]=   
    pT()*cos(phi()) +     z() *parent->showerVariables()[0];
  children[0]->showerVariables()[1]=   
    pT()*sin(phi()) +     z() *parent->showerVariables()[1];
  children[1]->showerVariables()[0]= 
    - pT()*cos(phi()) + (1.-z())*parent->showerVariables()[0];
  children[1]->showerVariables()[1]= 
    - pT()*sin(phi()) + (1.-z())*parent->showerVariables()[1];
  for(unsigned int ix=0;ix<2;++ix)
    children[ix]->showerVariables()[2]=
      sqrt(sqr(children[ix]->showerVariables()[0])+
	   sqr(children[ix]->showerVariables()[1]));
  // set up the colour connections
  splittingFn()->colourConnection(parent,children[0],children[1],partnerType,false);
  // make the products children of the parent
  parent->addChild(children[0]);
  parent->addChild(children[1]);
}

void FS_QTildeShowerKinematics1to2::
reconstructParent(const tShowerParticlePtr parent, 
	     const ParticleVector & children ) const {
  if(children.size() != 2) 
    throw Exception() << "FS_QTildeShowerKinematics1to2::updateParent() " 
		      << "Warning! too many children!" 
		      << Exception::eventerror;
  ShowerParticlePtr c1 = dynamic_ptr_cast<ShowerParticlePtr>(children[0]);
  ShowerParticlePtr c2 = dynamic_ptr_cast<ShowerParticlePtr>(children[1]);
  parent->showerParameters()[1]= 
    c1->showerParameters()[1] + c2->showerParameters()[1]; 
  parent->set5Momentum( c1->momentum() + c2->momentum() );
}

void FS_QTildeShowerKinematics1to2::reconstructLast(const tShowerParticlePtr theLast,
						    unsigned int iopt, 
						    Energy mass) const {
  // set beta component and consequently all missing data from that,
  // using the nominal (i.e. PDT) mass.
  Energy theMass = mass > ZERO  ?  mass : theLast->data().constituentMass();
  theLast->showerParameters()[1]=
    (sqr(theMass) + sqr(theLast->showerVariables()[2]) 
     - sqr( theLast->showerParameters()[0] )*pVector().m2())
    / ( 2.*theLast->showerParameters()[0]*p_dot_n() );
  // set that new momentum
  theLast->set5Momentum(sudakov2Momentum( theLast->showerParameters()[0],
					  theLast->showerParameters()[1], 
					  theLast->showerVariables()[0],
					  theLast->showerVariables()[1],iopt));
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
						 ShowerPartnerType::Type, 
						 bool) const {
  IdList ids(3);
  ids[0] = parent->id();
  ids[1] = children[0]->id();
  ids[2] = children[1]->id();
  vector<Energy> virtualMasses = SudakovFormFactor()->virtualMasses(ids);
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
  // set the values
  children[0]->showerVariables()[0]=   
    pT()*cos(phi()) +     z() *parent->showerVariables()[0];
  children[0]->showerVariables()[1]=   
    pT()*sin(phi()) +     z() *parent->showerVariables()[1];
  children[1]->showerVariables()[0]= 
    - pT()*cos(phi()) + (1.-z())*parent->showerVariables()[0];
  children[1]->showerVariables()[1]= 
    - pT()*sin(phi()) + (1.-z())*parent->showerVariables()[1];
  for(unsigned int ix=0;ix<2;++ix)
    children[ix]->showerVariables()[2]=
      sqrt(sqr(children[ix]->showerVariables()[0])+
	   sqr(children[ix]->showerVariables()[1]));
  for(unsigned int ix=0;ix<children.size();++ix) {
    if(children[ix]->children().empty()) continue;
    ShowerParticleVector newChildren;
    for(unsigned int iy=0;iy<children[ix]->children().size();++iy)
      newChildren.push_back(dynamic_ptr_cast<ShowerParticlePtr>
			    (children[ix]->children()[iy]));
    children[ix]->showerKinematics()->resetChildren(children[ix],newChildren);
  }
}
