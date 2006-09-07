// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbarAnalysis class.
//

#include "TTbarAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TTbarAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

TTbarAnalysis::~TTbarAnalysis() {}

namespace {
  inline Lorentz5Momentum getMomentum(tcPPtr particle) {
    return particle->momentum();
  } 
}


void TTbarAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {

  Lorentz5Momentum ptop, ptbar, ppair;  
  bool foundt = false;
  bool foundtbar = false;
  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));

  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if((**it).id() == ParticleID::t) {
      ptop = getMomentum(*it);
      foundt = true;
    } else if((**it).id() == ParticleID::tbar) {
      ptbar = getMomentum(*it);
      foundtbar = true;
    }
  }

  if (foundt && foundtbar) {
    ppair = ptop + ptbar;
    _pttop += ptop.perp()/GeV;
    _pttbar += ptbar.perp()/GeV;
    _ptpair += ppair.perp()/GeV;
    _ettop += ptop.et()/GeV;
    _ettbar += ptbar.et()/GeV;
    _etpair += ppair.et()/GeV;
    _etop += ptop.e()/GeV;
    _etbar += ptbar.e()/GeV;
    _epair += ppair.e()/GeV;
    _raptop += ptop.rapidity();
    _raptbar += ptbar.rapidity();
    _rappair += ppair.rapidity();
    _phitop += ptop.phi();
    _phitbar += ptbar.phi();
    _deltaphi += (ptop.vect()).deltaPhi(ptbar.vect());
    _mpair += ppair.m()/GeV;
    _etsum += (ptop.et() + ptbar.et())/GeV;
    _ptsum += (ptop.perp() + ptbar.perp())/GeV;
  } else {
    cerr << "Analysis/TTbarAnalysis: did not find ttbar pair in event " 
	 << event->number()  << ".\n";
    generator()->log() << "Analysis/TTbarAnalysis: " 
		       << "Found no ttbar pair in event " 
		       << event->number()  << ".\n"
		       << *event;    
  }  
}

LorentzRotation TTbarAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void TTbarAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
}

void TTbarAnalysis::analyze(tPPtr) {}

void TTbarAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void TTbarAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<TTbarAnalysis> TTbarAnalysis::initTTbarAnalysis;
// Definition of the static class description member.

void TTbarAnalysis::Init() {

  static ClassDocumentation<TTbarAnalysis> documentation
    ("There is no documentation for the TTbarAnalysis class");

}

