// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CKKWee2qq class.
//

#include "CKKWee2qq.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CKKWee2qq.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <cassert>

using namespace Herwig;

CKKWee2qq::~CKKWee2qq() {}

void CKKWee2qq::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _useColour;
}

void CKKWee2qq::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _useColour;
}

ClassDescription<CKKWee2qq> CKKWee2qq::initCKKWee2qq;
// Definition of the static class description member.

void CKKWee2qq::Init() {

  static ClassDocumentation<CKKWee2qq> documentation
    ("The e+e- to q qbar hard process.");


  static Switch<CKKWee2qq,bool> interfaceUseColour
    ("UseColour",
     "Use colour information for determining the hard process.",
     &CKKWee2qq::_useColour, true, false, false);
  static SwitchOption interfaceUseColourUseColourOn
    (interfaceUseColour,
     "On",
     "Use colour information.",
     true);
  static SwitchOption interfaceUseColourUseColourOff
    (interfaceUseColour,
     "Off",
     "Do not use colour information.",
     false);


}

bool CKKWee2qq::reachedHard (const vector<ClusteringParticleData>& particles) const {

  if (particles.size() != 4) return false;

  unsigned int goteMinus = 0;
  unsigned int gotePlus = 0;
  unsigned int gotQ = 0;
  unsigned int gotQbar = 0;

  unsigned int gotall = 0;

  for (unsigned int i = 0; i< 4; ++i) {
    if (particles[i].partonId.PDGId == 11 &&
	particles[i].partonId.state == ClusteringParticleState::initial) {
      goteMinus = i;
      gotall += 1;
    }
    if (particles[i].partonId.PDGId == -11 &&
	particles[i].partonId.state == ClusteringParticleState::initial) {
      gotePlus = i;
      gotall += 1;
    }
    if (abs(particles[i].partonId.PDGId) <7 && particles[i].partonId.PDGId > 0 &&
	particles[i].partonId.state == ClusteringParticleState::final) {
      gotQ = i;
      gotall += 1;
    }
    if (abs(particles[i].partonId.PDGId) <7 && particles[i].partonId.PDGId < 0 &&
	particles[i].partonId.state == ClusteringParticleState::final) {
      gotQbar = i;
      gotall += 1;
    }
  }

  if (gotall != 4) return false;

  if (particles[gotQ].partonId.PDGId + particles[gotQbar].partonId.PDGId != 0) return false;

  if (_useColour && particles[gotQ].colour != particles[gotQbar].antiColour) return false;

  return true;

}

