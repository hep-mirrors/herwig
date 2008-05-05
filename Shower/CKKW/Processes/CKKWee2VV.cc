// -*- C++ -*-
//
// CKKWee2VV.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CKKWee2VV class.
//

#include "CKKWee2VV.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CKKWee2VV.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <cassert>

using namespace Herwig;

CKKWee2VV::~CKKWee2VV() {}

void CKKWee2VV::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void CKKWee2VV::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<CKKWee2VV> CKKWee2VV::initCKKWee2VV;
// Definition of the static class description member.

void CKKWee2VV::Init() {

  static ClassDocumentation<CKKWee2VV> documentation
    ("The e+e- to V V hard process.");

}

bool CKKWee2VV::reachedHard (const vector<ClusteringParticleData>& particles) const {

  if (particles.size() != 4) return false;

  unsigned int goteMinus = 0;
  unsigned int gotePlus = 0;

  unsigned int gotAll = 0;

  for (unsigned int i = 0; i< 4; ++i) {
    if (particles[i].partonId.PDGId == 11 &&
	particles[i].partonId.state == ClusteringParticleState::initial) {
      goteMinus = i;
      gotAll +=1;
    }
    if (particles[i].partonId.PDGId == -11 &&
	particles[i].partonId.state == ClusteringParticleState::initial) {
      gotePlus = i;
      gotAll +=1;
    }

  }

  if (gotAll != 2) return false;

  long v1Id=0; long v2Id =0;

  for (unsigned int i = 0; i< 4; ++i)
    if (i != goteMinus || i != gotePlus) {
      if (v1Id ==0) v1Id = particles[i].partonId.PDGId;
      else v2Id = particles[i].partonId.PDGId;
    }
  if (particles[v1Id].partonId.state != ClusteringParticleState::final ||
      particles[v2Id].partonId.state != ClusteringParticleState::final) return false;

  if (abs(v1Id) == 24 && v1Id+v2Id ==0) return true;

  if (v1Id == 23 && v1Id == v2Id) return true;

  return false;

}


