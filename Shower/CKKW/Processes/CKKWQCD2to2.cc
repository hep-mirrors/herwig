// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CKKWQCD2to2 class.
//

#include "CKKWQCD2to2.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/UseRandom.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CKKWQCD2to2.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

CKKWQCD2to2::~CKKWQCD2to2() {}

void CKKWQCD2to2::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _useColour << _alpha;
}

void CKKWQCD2to2::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _useColour >> _alpha;
}

ClassDescription<CKKWQCD2to2> CKKWQCD2to2::initCKKWQCD2to2;
// Definition of the static class description member.

void CKKWQCD2to2::Init() {

  static ClassDocumentation<CKKWQCD2to2> documentation
    ("QCD 2 to 2 scatterings as hard proceswses in ME/PS merging.");

  static Switch<CKKWQCD2to2,bool> interfaceUseColour
    ("UseColour",
     "Use colour information for determining the hard process.",
     &CKKWQCD2to2::_useColour, true, false, false);
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

  
  static Reference<CKKWQCD2to2,ShowerAlpha> interfaceAlphaS
    ("AlphaS",
     "Alpha_s used for the hard process",
     &CKKWQCD2to2::_alpha, false, false, true, false, false);


}

bool CKKWQCD2to2::colourConservation (const vector<ClusteringParticleData>& particles) const {

  // check for colour conservation
  // by weighting the incoming colour with +1,
  // incoming anticolour with -1,
  // outgoing colour with -1,
  // outgoing anticolour with +1

  unsigned long coloursum = 0;

  coloursum += particles[0].colour - particles[0].antiColour;
  coloursum += particles[1].colour - particles[1].antiColour;

  coloursum -= particles[2].colour - particles[2].antiColour;
  coloursum -= particles[3].colour - particles[3].antiColour;

  return coloursum == 0;
  
}

bool CKKWQCD2to2::reachedHard (const vector<ClusteringParticleData>& particles) const {

  if (particles.size() != 4) return false;

  vector<ClusteringParticleData> sorted (4);

  bool gotFirstIn = false;
  bool gotFirstOut = false;

  for (vector<ClusteringParticleData>::const_iterator p = particles.begin(); 
       p != particles.end(); ++p) {
    if (p->partonId.state == ClusteringParticleState::initial) {
      if (gotFirstIn) {
	sorted[1] = *p;
      } else {
	gotFirstIn = true;
	sorted[0] = *p;
      }
    }
    if (p->partonId.state == ClusteringParticleState::final) {
      if (gotFirstOut) {
	sorted[3] = *p;
      } else {
	gotFirstOut = true;
	sorted[2] = *p;
      }
    }
  }

  // sort such that our processes read
  // q qbar -> X, g q -> X etc.

  if (sorted[0].partonId.PDGId < sorted[1].partonId.PDGId)
    swap(sorted[0],sorted[1]);
  if (sorted[2].partonId.PDGId < sorted[3].partonId.PDGId)
    swap(sorted[2],sorted[3]);

  // q qbar initiated

  if (abs(sorted[0].partonId.PDGId) < 7 && sorted[0].partonId.PDGId + sorted[1].partonId.PDGId == 0) {

    // q qbar final state

    if (abs(sorted[2].partonId.PDGId) < 7 && sorted[2].partonId.PDGId + sorted[3].partonId.PDGId == 0) {

      // q qbar -> q qbar

      if (abs(sorted[0].partonId.PDGId) == abs(sorted[2].partonId.PDGId)) {
	if (!_useColour) return true;
	else {

	  // s channel colour flow
	  if (sorted[0].colour == sorted[2].colour &&
	      sorted[1].antiColour == sorted[3].antiColour) return true;
	  
	  // t channel colour flow
	  if (sorted[0].colour == sorted[1].antiColour &&
	      sorted[2].colour == sorted[3].antiColour) return true;

	}

	return false;

      }

      // q qbar -> Q Qbar
      else {
	if (!_useColour) return true;
	else {

	  // only s channel flow here
	  if (sorted[0].colour == sorted[2].colour &&
	      sorted[1].antiColour == sorted[3].antiColour) return true;
	}

	return false;

      }

    }

    // gg final state

    if (sorted[2].partonId.PDGId == 21 && sorted[3].partonId.PDGId == 21) {
      if (!_useColour) return true;
      else

	// s channel and t channel flows are the same
      
	return sorted[0].colour != sorted[1].antiColour && colourConservation(sorted);

    }

  } // q qbar initiated


  // quark (anti)quark scattering

  if (sorted[0].partonId.PDGId == sorted[2].partonId.PDGId &&
      sorted[1].partonId.PDGId == sorted[3].partonId.PDGId &&
      abs(sorted[0].partonId.PDGId) < 7 &&
      abs(sorted[1].partonId.PDGId) < 7 &&
      sorted[0].partonId.PDGId != sorted[1].partonId.PDGId) {

    if (!_useColour) return true;
    else {

      if (sorted[0].partonId.PDGId * sorted[1].partonId.PDGId > 0) {

	// quark quark
	if (sorted[0].partonId.PDGId > 0 &&
	    sorted[0].colour == sorted[3].colour &&
	    sorted[1].colour == sorted[2].colour) return true;
	
	// antiquark antiquark
	if (sorted[0].partonId.PDGId < 0 &&
	    sorted[0].antiColour == sorted[3].antiColour &&
	    sorted[1].antiColour == sorted[2].antiColour) return true;
	
      } else {

	// quark antiquark
	if (sorted[0].colour == sorted[1].antiColour &&
	    sorted[2].colour == sorted[3].antiColour) return true;

      }

    }

    return false;

  } // quark (anti)quark scattering


  // g g initiated

  // gg -> gg
  if (sorted[0].partonId.PDGId == 21 &&
      sorted[1].partonId.PDGId == 21 &&
      sorted[2].partonId.PDGId == 21 &&
      sorted[3].partonId.PDGId == 21) {
    if (!_useColour) return true;
    else return colourConservation(sorted);
  }

  // gg -> q qbar
  if (sorted[0].partonId.PDGId == 21 &&
      sorted[1].partonId.PDGId == 21 &&
      abs(sorted[2].partonId.PDGId) < 7 &&
      sorted[2].partonId.PDGId + sorted[3].partonId.PDGId == 0) {
    if (!_useColour) return true;
    else

      // s channel and t channel flows are the same
      
      return sorted[2].colour != sorted[3].antiColour && colourConservation(sorted);
  }

  // gq -> gq
  if (sorted[0].partonId.PDGId == 21 && sorted[2].partonId.PDGId == 21 &&
      sorted[1].partonId.PDGId == sorted[3].partonId.PDGId &&
      abs(sorted[1].partonId.PDGId) < 7) {
    if (!_useColour) return true;
    else {
      if (sorted[1].partonId.PDGId > 0)
	return sorted[1].colour != sorted[3].colour && colourConservation(sorted);
      else
	return sorted[1].antiColour != sorted[3].antiColour && colourConservation(sorted);
    }
  }

  return false;

}

double CKKWQCD2to2::hardCouplings (const vector <tClusteringParticlePtr>& particles,double MEAlpha) {

  vector<tClusteringParticlePtr> sorted (4);

  bool gotFirstIn = false;
  bool gotFirstOut = false;

  for (vector<tClusteringParticlePtr>::const_iterator p = particles.begin(); 
       p != particles.end(); ++p) {
    if ((**p).pData().partonId.state == ClusteringParticleState::initial) {
      if (gotFirstIn) {
	sorted[1] = *p;
      } else {
	gotFirstIn = true;
	sorted[0] = *p;
      }
    }
    if ((**p).pData().partonId.state == ClusteringParticleState::final) {
      if (gotFirstOut) {
	sorted[3] = *p;
      } else {
	gotFirstOut = true;
	sorted[2] = *p;
      }
    }
  }

  // get the mandelstam invariants

  Lorentz5Momentum ps = sorted[0]->momentum()+sorted[1]->momentum();
  Lorentz5Momentum pt = sorted[0]->momentum()-sorted[2]->momentum();
  Lorentz5Momentum pu = sorted[0]->momentum()-sorted[3]->momentum();

  Energy2 s = ps*ps;
  Energy2 t = pt*pt;
  Energy2 u = pu*pu;

  return sqr(_alpha->value(2.*s*t*u/(sqr(s)+sqr(t)+sqr(u)))/MEAlpha);

}
