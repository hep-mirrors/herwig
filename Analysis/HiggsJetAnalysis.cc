// -*- C++ -*-
//
// HiggsJetAnalysis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiggsJetAnalysis class.
//

#include "HiggsJetAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HiggsJetAnalysis::HiggsJetAnalysis() :
  _pth(0.,250.,100), _pthZoom(35.,65.,100), 
  _raph(-10.,10.,100), _phih(-Constants::pi,Constants::pi,100) {}

void HiggsJetAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  _pth.topdrawOutput(outfile,Frame,"BLACK","pt of Higgs");
  _pth.topdrawOutput(outfile,Frame|Ylog,"BLACK","pt of Higgs");
  _pthZoom.topdrawOutput(outfile,Frame,"BLACK","35<pt/GeV<65 of Higgs");
  _pthZoom.topdrawOutput(outfile,Frame|Ylog,"BLACK","35<pt/GeV<65 of Higgs");
  _raph.topdrawOutput(outfile,Frame,"BLACK","Rapidity of h");
  _raph.topdrawOutput(outfile,Frame|Ylog,"BLACK","Rapidity of h");
  _phih.topdrawOutput(outfile,Frame,"BLACK","Azimuth of h");
  _phih.topdrawOutput(outfile,Frame|Ylog,"BLACK","Azimuth of h");
  outfile.close();
}

namespace {
  bool isLastInShower(const Particle & p) {
    return p.children().size() > 1 
      && p.children()[0]->id() != p.id()
      && p.children()[1]->id() != p.id();
  }

  struct Higgs {
    static bool AllCollisions() { return false; }
    static bool AllSteps() { return true; }
    // ===
    // pick the last instance from the shower
    static bool FinalState() { return false; }
    static bool Intermediate() { return true; }
    // ===
    static bool Check(const Particle & p) { 
      return p.id() == ParticleID::h0 && isLastInShower(p);
    }
  };


}


void HiggsJetAnalysis::analyze(tEventPtr event, long, int, int) {
  tcParticleSet higgses;
  event->select(inserter(higgses), ThePEG::ParticleSelector<Higgs>());

  if ( higgses.empty() )
    return;
  else if ( higgses.size() > 1 ) {
    cerr << "\nMultiple h0 found. Only binning first one.\n";
  }

  tcPPtr higgs = *higgses.begin();

  Lorentz5Momentum ph = higgs->momentum();
  double pt = ph.perp()/GeV;
  (_pth)+=(pt);
  (_pthZoom)+=(pt);
  double rap = 0.5*log((ph.e()+ph.z())/(ph.e()-ph.z()));
  (_raph)+=(rap);
  (_phih)+=ph.phi();
}

NoPIOClassDescription<HiggsJetAnalysis> HiggsJetAnalysis::initHiggsJetAnalysis;
// Definition of the static class description member.

void HiggsJetAnalysis::Init() {

  static ClassDocumentation<HiggsJetAnalysis> documentation
    ("Standard analysis of a single h0 after showering.");

}

