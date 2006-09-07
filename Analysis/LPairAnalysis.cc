// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LPairAnalysis class.
//

#include "LPairAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LPairAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LPairAnalysis::~LPairAnalysis() {}

namespace {
  inline Lorentz5Momentum getMomentum(tcPPtr particle) {
    return particle->momentum();
  }

  inline bool isLeptonPlus(tcPPtr p) {
    if ( p->id() == ParticleID::eplus ||
	 p->id() == ParticleID::muplus ||
	 p->id() == ParticleID::tauplus ) {
      return true;
    }else {
      return false;
    }
  }

  inline bool isLeptonMinus(tcPPtr p) {
    if ( p->id() == ParticleID::eminus ||
	 p->id() == ParticleID::muminus ||
	 p->id() == ParticleID::tauminus ) {
      return true;
    }else {
      return false;
    }
  }
  inline bool isFromTop(tcPPtr p) {
    while (p->parents()[0] && p->parents().size() == 1) {
      p = p->parents()[0];
      if (abs(p->id()) == ParticleID::t) {
	return true;
      }
    }
    return false; 
  } 
}


void LPairAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  Lorentz5Momentum ppair, plp, plm;  
  bool foundlp = false;
  bool foundlm = false;
  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));

  // find highest pt lepton+ and lepton- in the event resp.
  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if( isLeptonPlus(*it) ) {
      //     if ( getMomentum(*it).perp() > plp.perp() ) {
      if (isFromTop(*it)) {
	plp = getMomentum(*it);
	foundlp = true;
      }
    } else if( isLeptonMinus(*it) ) {
      //      if ( getMomentum(*it).perp() > plm.perp() ) {
      if (isFromTop(*it)) {
	plm = getMomentum(*it);
	foundlm = true;
      }
    }
  }
  
  if (foundlp && foundlm) {
    ppair = plp + plm;
    _ptp += plp.perp()/GeV;
    _ptm += plm.perp()/GeV;
    _ptpair += ppair.perp()/GeV;
    _etp += plp.et()/GeV;
    _etm += plm.et()/GeV;
    _etpair += ppair.et()/GeV;
    _ep += plp.e()/GeV;
    _em += plm.e()/GeV;
    _epair += ppair.e()/GeV;
    _rapp += plp.rapidity();
    _rapm += plm.rapidity();
    _rappair += ppair.rapidity();
    _phip += plp.phi();
    _phim += plm.phi();
    _deltaphi += (plp.vect()).deltaPhi(plm.vect());
    _mpair += ppair.m()/GeV;
    _etsum += (plp.et() + plm.et())/GeV;
    _ptsum += (plp.perp() + plm.perp())/GeV;
  } else {
    cerr << "Analysis/LPairAnalysis: did not find suitable lepton"
	 << " pair in event " << event->number()  << " ("
	 << (foundlp ? "+" : "0") 
	 << (foundlm ? "-" : "0") 
	 << ").\n";
    generator()->log() << "Analysis/LPairAnalysis: " 
		       << "Found no suitable lepton pair in event " 
		       << event->number()  << ".\n"
		       << *event;    
  }  
}

LorentzRotation LPairAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void LPairAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
}

void LPairAnalysis::analyze(tPPtr) {}

void LPairAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void LPairAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<LPairAnalysis> LPairAnalysis::initLPairAnalysis;
// Definition of the static class description member.

void LPairAnalysis::Init() {

  static ClassDocumentation<LPairAnalysis> documentation
    ("There is no documentation for the LPairAnalysis class");

}

