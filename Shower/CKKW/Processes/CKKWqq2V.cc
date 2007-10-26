// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CKKWqq2V class.
//

#include "CKKWqq2V.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CKKWqq2V.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <cassert>

using namespace Herwig;

CKKWqq2V::~CKKWqq2V() {}

void CKKWqq2V::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _useColour;
}

void CKKWqq2V::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _useColour;
}

ClassDescription<CKKWqq2V> CKKWqq2V::initCKKWqq2V;
// Definition of the static class description member.

void CKKWqq2V::Init() {

  static ClassDocumentation<CKKWqq2V> documentation
    ("The q q to V hard process.");


  static Switch<CKKWqq2V,bool> interfaceUseColour
    ("UseColour",
     "Use colour information for determining the hard process.",
     &CKKWqq2V::_useColour, true, false, false);
  static SwitchOption interfaceUseColourUseColourOn
    (interfaceUseColour,
     "Yes",
     "Use colour information.",
     true);
  static SwitchOption interfaceUseColourUseColourOff
    (interfaceUseColour,
     "No",
     "Do not use colour information.",
     false);


}

bool CKKWqq2V::reachedHard (const vector<ClusteringParticleData>& particles) const {

  if (particles.size() != 3) return false;

  unsigned int gotQ = 0;
  unsigned int gotQbar = 0;
  unsigned int gotV = 0;

  unsigned int gotAll = 0;

  for (unsigned int i = 0; i<3 ; ++i) {
    if ((particles[i].partonId.PDGId == 23 || abs(particles[i].partonId.PDGId) == 24) &&
	particles[i].partonId.state == ClusteringParticleState::final) {
      gotV = i;
      gotAll += 1;
    }
    if (abs(particles[i].partonId.PDGId) <7 && particles[i].partonId.PDGId > 0 &&
	particles[i].partonId.state == ClusteringParticleState::initial) {
      gotQ = i;
      gotAll += 1;
    }
    if (abs(particles[i].partonId.PDGId) <7 && particles[i].partonId.PDGId < 0 &&
	particles[i].partonId.state == ClusteringParticleState::initial) {
      gotQbar = i;
      gotAll += 1;
    }
  }

  if (gotAll != 3) return false;

  if (_useColour && particles[gotQ].colour != particles[gotQbar].antiColour) return false;

  if(particles[gotV].partonId.PDGId == 23) {
    if (particles[gotQ].partonId.PDGId + particles[gotQbar].partonId.PDGId != 0) return false;
    else return true;
  }

  if(particles[gotV].partonId.PDGId == 24) {
    if (particles[gotQ].partonId.PDGId + particles[gotQbar].partonId.PDGId != 1) return false;
    else return true;
  }

  if(particles[gotV].partonId.PDGId == -24) {
    if (particles[gotQ].partonId.PDGId + particles[gotQbar].partonId.PDGId != -1) return false;
    else return true;
  }

  return false;

}


