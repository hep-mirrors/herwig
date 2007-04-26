// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpinCorrAnalysis class.
//

#include "SpinCorrAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;

SpinCorrAnalysis::~SpinCorrAnalysis() {}

void SpinCorrAnalysis::analyze(tEventPtr event, long, int, int) {
  //  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  ParticleVector pzero = event->primarySubProcess()->outgoing();
  // find the left squark
  PPtr bsmquark;
  for(ParticleVector::const_iterator it = pzero.begin(); it != pzero.end(); 
      ++it) {
    if((**it).id() == 1000002 || (**it).id() == 5100002 ) {
      bsmquark = *it;
      break;
    }
  }
  if(!bsmquark) return;
  tPPtr quark, lepton, kk1_lept;
  if(bsmquark->id() == 1000002) {
    tPPtr  chi2 = findChild(bsmquark, 1000023);
    quark = findChild(bsmquark, 2);
    
    if(!quark || !chi2) return;

    lepton = findChild(chi2, 11);
    tPPtr slepton = findChild(chi2, 2000011);
  
    if(!slepton || !lepton) return;
  }
  else {
    quark = findChild(bsmquark, 2);
    tPPtr kkz1 = findChild(bsmquark, 5100023);

    if(!quark || !kkz1) return;
    lepton = findChild(kkz1, 11);
    kk1_lept = findChild(kkz1, 6100011);

    if(!lepton || !kk1_lept) return;
      
  }
  if(lepton && quark) {
    Lorentz5Momentum ptotal = quark->momentum() + lepton->momentum();
    Energy invmass = sqrt(ptotal.m2());
    
    if(lepton->id() < 0) {
      *thelPlus += invmass/theMaxInvMass;
    }
    else {
      *thelMinus += invmass/theMaxInvMass;
    }
}
    
}

LorentzRotation SpinCorrAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
}

void SpinCorrAnalysis::analyze(const tPVector &) {
}

void SpinCorrAnalysis::analyze(tPPtr) {}

void SpinCorrAnalysis::persistentOutput(PersistentOStream & os) const {
  os << theMaxInvMass;
}

void SpinCorrAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> theMaxInvMass;
}

ClassDescription<SpinCorrAnalysis> SpinCorrAnalysis::initSpinCorrAnalysis;
// Definition of the static class description member.

void SpinCorrAnalysis::Init() {

  static ClassDocumentation<SpinCorrAnalysis> documentation
    ("There is no documentation for the SpinCorrAnalysis class");

}

tPPtr SpinCorrAnalysis::findDecayingParticle(tPPtr parent) const {
  ParticleVector children = parent->children();
  ParticleVector::size_type nprod = children.size();
  //  if(nprod == 0) return parent;
  bool final(true);
  tPPtr copy;
  for(unsigned int i = 0; i < nprod; ++i) {
    if(children[i]->id() == parent->id()) {
      final = false;
      copy = children[i];
      break;
    }
  }
  if(!final) parent = findDecayingParticle(copy);
  return parent;
}


tPPtr SpinCorrAnalysis::findChild(tPPtr parent, long child) const {
  //cerr << "testing A \n";
  //find children after shower
  if(!parent) return tPPtr();
  //  parent =  findDecayingParticle(parent);
  ParticleVector children = parent->children();
  tPPtr childparticle;
  for(unsigned int i = 0; i < children.size(); ++i) {
    //    cerr << "children: " << children[i]->PDGName() << endl;
    //cerr << i << "\n";
    if(abs(children[i]->id()) == child) {
      childparticle = children[i];
      break;
    }
  }

  //  cerr << "testing D \n";
  return childparticle;
}
