// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QCDEmitterSpectatorClusterer class.
//

#include "QCDEmitterSpectatorClusterer.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QCDEmitterSpectatorClusterer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

QCDEmitterSpectatorClusterer::~QCDEmitterSpectatorClusterer() {}

void QCDEmitterSpectatorClusterer::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void QCDEmitterSpectatorClusterer::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

AbstractClassDescription<QCDEmitterSpectatorClusterer> QCDEmitterSpectatorClusterer::initQCDEmitterSpectatorClusterer;
// Definition of the static class description member.

void QCDEmitterSpectatorClusterer::Init() {

  static ClassDocumentation<QCDEmitterSpectatorClusterer> documentation
    ("QCD emitter spectator clusterings.");

}

vector<ClusteringConfigurationPtr>
QCDEmitterSpectatorClusterer::configurations (const vector<ClusteringParticleData>& in) {
  vector<ClusteringParticleData> out;
  vector<ClusteringConfigurationPtr> result;

  ClusteringParticleData em;

  bool isQCD = false;

  // (0 1) 2

  em = emergingLine(in[0],in[1],isQCD);
  if (isQCD && colourConnected(em,in[2])) {
    out.push_back(em);
    out.push_back(in[2]);
    result.push_back(new_ptr(EmitterSpectatorConfiguration(in,2,out,1,ClusteringInteractionType::QCD,this)));
  }
  out.clear();

  // (1 2) 0

  em = emergingLine(in[1],in[2],isQCD);
  if (isQCD && colourConnected(em,in[0])) {
    out.push_back(em);
    out.push_back(in[0]);
    result.push_back(new_ptr(EmitterSpectatorConfiguration(in,0,out,1,ClusteringInteractionType::QCD,this)));
  }
  out.clear();

  // (0 2) 1

  em = emergingLine(in[0],in[2],isQCD);
  if (isQCD && colourConnected(em,in[1])) {
    out.push_back(em);
    out.push_back(in[1]);
    result.push_back(new_ptr(EmitterSpectatorConfiguration(in,1,out,1,ClusteringInteractionType::QCD,this)));
  }
  out.clear();

  return result;

}
