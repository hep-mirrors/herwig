// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EmitterSpectatorConfiguration class.
//

#include "EmitterSpectatorConfiguration.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EmitterSpectatorConfiguration.tcc"
#endif


using namespace Herwig;

EmitterSpectatorConfiguration::~EmitterSpectatorConfiguration() {}

NoPIOClassDescription<EmitterSpectatorConfiguration> EmitterSpectatorConfiguration::initEmitterSpectatorConfiguration;
// Definition of the static class description member.

void EmitterSpectatorConfiguration::Init() {

  static ClassDocumentation<EmitterSpectatorConfiguration> documentation
    ("EmitterSpectatorConfiguration stores quantum number assignments for "
     "clusterings in an emitter-spectator scenario.");

}

#ifdef HERWIG_DEBUG_CKKW

void EmitterSpectatorConfiguration::debugDump (ostream& os) {
  ClusteringConfiguration::debugDump(os);
  os << "spectator before clustering " << _childSpectator
     << " after " << _parentSpectator << endl;
}

#endif 
