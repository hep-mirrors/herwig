// -*- C++ -*-
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

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HiggsJetAnalysis.tcc"
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

LorentzRotation HiggsJetAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void HiggsJetAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
}

void HiggsJetAnalysis::analyze(tPPtr) {}

NoPIOClassDescription<HiggsJetAnalysis> HiggsJetAnalysis::initHiggsJetAnalysis;
// Definition of the static class description member.

void HiggsJetAnalysis::Init() {

  static ClassDocumentation<HiggsJetAnalysis> documentation
    ("Standard analysis of a single h0 after showering.");

}

