// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EWKEmitterSpectatorClusterer class.
//

#include "EWKEmitterSpectatorClusterer.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EWKEmitterSpectatorClusterer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

EWKEmitterSpectatorClusterer::~EWKEmitterSpectatorClusterer() {}

void EWKEmitterSpectatorClusterer::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void EWKEmitterSpectatorClusterer::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

AbstractClassDescription<EWKEmitterSpectatorClusterer> EWKEmitterSpectatorClusterer::initEWKEmitterSpectatorClusterer;
// Definition of the static class description member.

void EWKEmitterSpectatorClusterer::Init() {

  static ClassDocumentation<EWKEmitterSpectatorClusterer> documentation
    ("EWK emitter spectator clusterings.");

}

vector<ClusteringConfigurationPtr>
EWKEmitterSpectatorClusterer::configurations (const vector<ClusteringParticleData>& in) {
  vector<ClusteringParticleData> out;
  vector<ClusteringConfigurationPtr> result;

  ClusteringParticleData em;

  bool isEWK = false;

  // (0 1) 2

  em = emergingLine(in[0],in[1],isEWK);
  if (isEWK) {
    out.push_back(em);
    out.push_back(in[2]);
    result.push_back(new_ptr(EmitterSpectatorConfiguration(in,2,out,1,ClusteringInteractionType::EWK,this)));
  }
  out.clear();

  // (1 2) 0

  em = emergingLine(in[1],in[2],isEWK);
  if (isEWK) {
    out.push_back(em);
    out.push_back(in[0]);
    result.push_back(new_ptr(EmitterSpectatorConfiguration(in,0,out,1,ClusteringInteractionType::EWK,this)));
  }
  out.clear();

  // (0 2) 1

  em = emergingLine(in[0],in[2],isEWK);
  if (isEWK) {
    out.push_back(em);
    out.push_back(in[1]);
    result.push_back(new_ptr(EmitterSpectatorConfiguration(in,1,out,1,ClusteringInteractionType::EWK,this)));
  }
  out.clear();

  return result;

}
