// -*- C++ -*-
//
// EmitterSpectatorClustering.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EmitterSpectatorClustering class.
//

#include "EmitterSpectatorClustering.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/PDT/ParticleData.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EmitterSpectatorClustering.tcc"
#endif


using namespace Herwig;

EmitterSpectatorClustering::~EmitterSpectatorClustering() {}

NoPIOClassDescription<EmitterSpectatorClustering> EmitterSpectatorClustering::initEmitterSpectatorClustering;
// Definition of the static class description member.

void EmitterSpectatorClustering::Init() {

  static ClassDocumentation<EmitterSpectatorClustering> documentation
    ("3 -> 2 clustering involving an emitter-spectator scenario.");

}

void EmitterSpectatorClustering::generateSudakovBasis () {

  Lorentz5Momentum n = spectatorBeforeClustering()->momentum();
  n.setMass(0.*GeV); n.rescaleEnergy();

  Lorentz5Momentum p;

  if (emission().first->pData().partonId.state == ClusteringParticleState::final &&
      emission().second->pData().partonId.state == ClusteringParticleState::final) {
    p = emission().first->momentum() + emission().second->momentum();
  }

  if (emission().first->pData().partonId.state == ClusteringParticleState::initial &&
      emission().second->pData().partonId.state == ClusteringParticleState::final) {
    p = emission().first->momentum() - emission().second->momentum();
  }

  if (emission().first->pData().partonId.state == ClusteringParticleState::final &&
      emission().second->pData().partonId.state == ClusteringParticleState::initial) {
    p = - emission().first->momentum() + emission().second->momentum();
  }

  Energy emmMass = getParticleData(emitter()->pData().partonId.PDGId)->mass();

  p.setMass(emmMass); p.rescaleEnergy();

  _sudakovBasis = make_pair(p,n);

}

#ifdef HERWIG_DEBUG_CKKW

void EmitterSpectatorClustering::debugDump (ostream& os) {
  Clustering::debugDump(os);
  os << "spectator before clustering " << _childSpectator
     << " after " << _parentSpectator << endl;
}

#endif 
