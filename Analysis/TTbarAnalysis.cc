// -*- C++ -*-
//
// TTbarAnalysis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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

TTbarAnalysis::TTbarAnalysis() :
  _pttop(0.,350.,100), _pttbar(0.,350.,100), _ptpair(0.,350.,100), 
  _ettop(0.,350.,100), _ettbar(0.,350.,100), _etpair(0.,350.,100), 
  _etop(0.,3000.,100), _etbar(0.,3000.,100), _epair(0.,6000.,100), 
  _raptop(-5.,5.,100), _raptbar(-5.,5.,100), _rappair(-5.,5.,100), 
  _phitop  (-Constants::pi,Constants::pi,50), 
  _phitbar (-Constants::pi,Constants::pi,50), 
  _deltaphi(-Constants::pi,Constants::pi,100), _mpair(300,1500,100), 
  _etsum(0.,700.,100), _ptsum(0.,700.,100) {}

void TTbarAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  _pttop.topdrawOutput(outfile,Frame|Ylog,"RED","pt t, tbar");
  _pttbar.topdrawOutput(outfile,Ylog,"BLUE","pt tbar");
  _ptpair.topdrawOutput(outfile,Frame|Ylog,"BLACK","pt pair");
  _ettop.topdrawOutput(outfile,Frame|Ylog,"RED","Et t, tbar");
  _ettbar.topdrawOutput(outfile,Ylog,"BLUE","Et tbar");
  _etpair.topdrawOutput(outfile,Frame|Ylog,"BLACK","Et pair");
  _etop.topdrawOutput(outfile,Frame|Ylog,"RED","E t, tbar");
  _etbar.topdrawOutput(outfile,Ylog,"BLUE","E tbar");
  _epair.topdrawOutput(outfile,Frame|Ylog,"BLACK","E pair");
  _raptop.topdrawOutput(outfile,Frame,"RED","y t, tbar");
  _raptbar.topdrawOutput(outfile,None,"BLUE","y tbar");
  _rappair.topdrawOutput(outfile,Frame,"BLACK","y pair");
  _phitop.topdrawOutput(outfile,Frame,"RED","phi t, tbar");
  _phitbar.topdrawOutput(outfile,None,"BLUE","phi tbar");
  _deltaphi.topdrawOutput(outfile,Frame,"BLACK","Delta phi");
  _mpair.topdrawOutput(outfile,Frame|Ylog,"BLACK","M pair");
  _etsum.topdrawOutput(outfile,Frame|Ylog,"BLACK","pt sum");
  _ptsum.topdrawOutput(outfile,Frame|Ylog,"BLACK","Et sum");
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

NoPIOClassDescription<TTbarAnalysis> TTbarAnalysis::initTTbarAnalysis;
// Definition of the static class description member.

void TTbarAnalysis::Init() {

  static ClassDocumentation<TTbarAnalysis> documentation
    ("Standard analysis of a t/tbar pair after showering.");

}

