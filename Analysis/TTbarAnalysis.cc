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

namespace {

  bool isLastInShower(const Particle & p) {
    return p.children().size() > 1 
      && p.children()[0]->id() != p.id()
      && p.children()[1]->id() != p.id();
  }

  struct TTBar {
    static bool AllCollisions() { return false; }
    static bool AllSteps() { return true; }
    // ===
    // pick the last instance from the shower
    static bool FinalState() { return false; }
    static bool Intermediate() { return true; }
    // ===
    static bool Check(const Particle & p) { 
      return abs(p.id()) == ParticleID::t && isLastInShower(p);
    }
  };

}


void TTbarAnalysis::analyze(tEventPtr event, long, int, int) {

  Lorentz5Momentum ptop, ptbar, ppair;  
  bool foundt = false;
  bool foundtbar = false;

  tcParticleSet particles;
  event->select(inserter(particles), ThePEG::ParticleSelector<TTBar>());

  if ( particles.empty() )
    return;

  for(tcParticleSet::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if((**it).id() == ParticleID::t) {
      ptop = (*it)->momentum();
      foundt = true;
    } else if((**it).id() == ParticleID::tbar) {
      ptbar = (*it)->momentum();
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

LorentzRotation TTbarAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void TTbarAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
}

void TTbarAnalysis::analyze(tPPtr) {}

NoPIOClassDescription<TTbarAnalysis> TTbarAnalysis::initTTbarAnalysis;
// Definition of the static class description member.

void TTbarAnalysis::Init() {

  static ClassDocumentation<TTbarAnalysis> documentation
    ("Standard analysis of a t/tbar pair after showering.");

}

