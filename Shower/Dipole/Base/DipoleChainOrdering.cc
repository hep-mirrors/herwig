// -*- C++ -*-
//
// DipoleChainOrdering.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleChainOrdering class.
//

#include "DipoleChainOrdering.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DipoleChainOrdering::DipoleChainOrdering() 
  : DipoleEvolutionOrdering(), virtualityOrdering(false) {}

DipoleChainOrdering::~DipoleChainOrdering() {}

IBPtr DipoleChainOrdering::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleChainOrdering::fullclone() const {
  return new_ptr(*this);
}

Energy DipoleChainOrdering::hardScale(tPPtr emitter, tPPtr spectator,
				      double emitterX, double spectatorX,
				      const DipoleSplittingKernel& split,
				      const DipoleIndex& index) const {

  Energy scale = 
    split.splittingKinematics()->dipoleScale(emitter->momentum(),
					     spectator->momentum());

  return 
    virtualityOrdering ?
    split.splittingKinematics()->QMax(scale,emitterX,spectatorX,index,split) :
    split.splittingKinematics()->ptMax(scale,emitterX,spectatorX,index,split);

}


void DipoleChainOrdering::setEvolutionScale(Energy scale,
					    const DipoleSplittingInfo&,
					    DipoleChain& chain,
					    pair<list<Dipole>::iterator,list<Dipole>::iterator>) const {

  for ( list<Dipole>::iterator dip = chain.dipoles().begin();
	dip != chain.dipoles().end(); ++dip ) {

    if ( dip->emitterScale(make_pair(true,false)) > scale )
      dip->emitterScale(make_pair(true,false),scale);

    if ( dip->emitterScale(make_pair(false,true)) > scale )
      dip->emitterScale(make_pair(false,true),scale);

  }

}

void DipoleChainOrdering::setEvolutionScale(Energy scale,
					    const DipoleSplittingInfo&,
					    DipoleChain& chain,
					    list<Dipole>::iterator) const {

  for ( list<Dipole>::iterator dip = chain.dipoles().begin();
	dip != chain.dipoles().end(); ++dip ) {

    if ( dip->emitterScale(make_pair(true,false)) > scale )
      dip->emitterScale(make_pair(true,false),scale);

    if ( dip->emitterScale(make_pair(false,true)) > scale )
      dip->emitterScale(make_pair(false,true),scale);

  }

}

Energy DipoleChainOrdering::evolutionScale(const DipoleSplittingInfo& split,
					   const DipoleSplittingKernel& spkernel) const {
  return 
    virtualityOrdering ?
    spkernel.splittingKinematics()->QFromPt(split.lastPt(),split) :
    split.lastPt();
}

Energy DipoleChainOrdering::maxPt(Energy scale,
				  const DipoleSplittingInfo& split,
				  const DipoleSplittingKernel& spkernel) const {
  return 
    virtualityOrdering ?
    spkernel.splittingKinematics()->ptMax(scale,split.emitterX(),split.spectatorX(),split.index(),spkernel) :
    scale;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleChainOrdering::persistentOutput(PersistentOStream & os) const {
  os << virtualityOrdering;
}

void DipoleChainOrdering::persistentInput(PersistentIStream & is, int) {
  is >> virtualityOrdering;
}

ClassDescription<DipoleChainOrdering> DipoleChainOrdering::initDipoleChainOrdering;
// Definition of the static class description member.

void DipoleChainOrdering::Init() {

  static ClassDocumentation<DipoleChainOrdering> documentation
    ("DipoleChainOrdering performs ordering on "
     "complete colour singlet dipole chains.");


  static Switch<DipoleChainOrdering,bool> interfaceVirtualityOrdering
    ("Ordering",
     "[experimental] Switch between virtuality and pt ordering.",
     &DipoleChainOrdering::virtualityOrdering, false, false, false);
  static SwitchOption interfaceVirtualityOrderingPt
    (interfaceVirtualityOrdering,
     "Pt",
     "Perform pt ordering",
     false);
  static SwitchOption interfaceVirtualityOrderingVirtuality
    (interfaceVirtualityOrdering,
     "Virtuality",
     "Perform virtuality ordering",
     true);

  interfaceVirtualityOrdering.rank(-1);

}

