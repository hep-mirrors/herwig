// -*- C++ -*-
//
// This is the implementation of the inlined, non-templated member
// functions of the Cutter class.
//

#include "Cutter.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


#include "fastjet/ClusterSequence.hh"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Handlers/EventHandler.h"
#include<iostream>



using namespace arnold;

Cutter::Cutter()
  : StepHandler(), 
    theRparam(0.7),theEtaDetector(5.0),theJetDefPt(20) {}

Cutter::~Cutter() {}

inline double deltaphi(double phi1, double phi2){
  double diff=phi1-phi2;
  if(diff<-Constants::pi){
    diff+=(2.0*Constants::pi);
  }
  else if (diff>Constants::pi){
    diff-=(2.0*Constants::pi);
    }
  return diff;
}

bool Cutter::allowedParticle (const Particle &p) {
  if (p.id()!=12 && p.id()!=14 && p.id()!=16 && p.id()!=18 &&p.id()!=25 && p.id()!=36 && p.id()!=82) return true;
  return false;
}

vector<fastjet::PseudoJet> Cutter::recombinables(const ParticleVector& p, bool sort_out){
  vector<fastjet::PseudoJet> recombinables;
  for (ParticleVector::const_iterator iter=p.begin(); iter!=p.end(); iter++){
    if (sort_out && !allowedParticle(**iter))
      continue;
    recombinables.push_back(fastjet::PseudoJet( (**iter).momentum().x()/GeV,(**iter).momentum().y()/GeV,(**iter).momentum().z()/GeV, (**iter).momentum().t()/GeV ));
  }
  return recombinables;
}

vector<fastjet::PseudoJet> Cutter::recombinables(const tPVector& p, bool sort_out){
  vector<fastjet::PseudoJet> recombinables;
  for (tPVector::const_iterator iter=p.begin(); iter!=p.end(); iter++){
    if (sort_out && !allowedParticle(**iter))
      continue;
    recombinables.push_back(fastjet::PseudoJet( (**iter).momentum().x()/GeV,(**iter).momentum().y()/GeV,(**iter).momentum().z()/GeV, (**iter).momentum().t()/GeV ));
  }
  return recombinables;
}

vector<fastjet::PseudoJet> Cutter::recombine(const vector<fastjet::PseudoJet>& p){
  //use fastjet to recombine jets
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
  fastjet::JetDefinition jet_def(fastjet::kt_algorithm, theRparam, recomb_scheme, strategy);
  fastjet::ClusterSequence clust_seq(p, jet_def);
  
  vector<fastjet::PseudoJet> jets = clust_seq.inclusive_jets(theJetDefPt);
  jets = sorted_by_pt(jets);

  return jets;
}

vector<fastjet::PseudoJet> Cutter::getInRange(const vector<fastjet::PseudoJet>& j){
  vector<fastjet::PseudoJet> jets;
  for (vector<fastjet::PseudoJet>::const_iterator iter = j.begin();
       iter != j.end(); iter++){
    if ( abs(iter->pseudorapidity()) < theEtaDetector ) jets.push_back(*iter);
  }
  return jets;
}



void Cutter::
handle(EventHandler & eh, const tPVector & tagged,
       const Hint &) throw(Veto, Stop, Exception) {
  // Implement the Handle method here.
  // Note that if the method actually does anything to the current event
  // the changes should be inserted in a new step which should be obtained
  // by 'ch.newStep()'.
  // Note also that the general advice is to only consider the particles in
  // the 'tagged' vector.

  if (theNLOAnalysis->runsPlainNLO()) return;

  const ParticleVector& hard_partons = eh.currentEvent()->primarySubProcess()->outgoing();
  vector<fastjet::PseudoJet> partons = sorted_by_pt(recombinables(hard_partons,false));


  vector<fastjet::PseudoJet> input_particles = recombinables(tagged);
  vector<fastjet::PseudoJet> jets = getInRange(recombine(input_particles));
  vector<fastjet::PseudoJet> recombined = recombine(input_particles);

  if ( !theNLOAnalysis->passCuts(jets) ){
    if (theExcessEventsAnalysis) {
      theExcessEventsAnalysis->setJetCache(jets);
      theExcessEventsAnalysis->analyze(eh.currentEvent(),1,1,0);
    }
    throw Veto();
  }
  
  theNLOAnalysis->setJetCache(jets);

}

IBPtr Cutter::clone() const {
  return new_ptr(*this);
}

IBPtr Cutter::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void Cutter::persistentOutput(PersistentOStream & os) const {
  os << theNLOAnalysis<< theRparam << theEtaDetector
     << theJetDefPt << theExcessEventsAnalysis; 
}

void Cutter::persistentInput(PersistentIStream & is, int) {
  is >> theNLOAnalysis>> theRparam  >> theEtaDetector 
     >> theJetDefPt >> theExcessEventsAnalysis;
}

ClassDescription<Cutter> Cutter::initCutter;
// Definition of the static class description member.

void Cutter::Init() {

  typedef bool (arnold::Cutter::*IGFN)() const;
  typedef void (arnold::Cutter::*ISFN)(bool);
  typedef Energy (arnold::Cutter::*IGFNK)() const;
  typedef void (arnold::Cutter::*ISFNK)(Energy);

 

  static ClassDocumentation<Cutter> documentation
    ("There is no documentation for the Cutter class yet");

  static Reference<Cutter,NLOAnalysis> interfaceNLOAnalysis
    ("NLOAnalysis",
     "A reference to NLOAnalysis",
     &Cutter::theNLOAnalysis, true, false, true, false, false);

  static Reference<Cutter,NLOAnalysis> interfaceExcessEventsAnalysis
    ("ExcessEventsAnalysis",
     "A reference to ExcessEventsAnalysis",
     &Cutter::theExcessEventsAnalysis, true, false, true, true, false);

  static Parameter<Cutter,double> interfaceRparam
    ("Rparam",
     "The R parameter of the jet definition algorithm",
     &Cutter::theRparam, 0.7, 0.0, 2.0, true, false, true);

  static Parameter<Cutter,double> interfacey_max
    ("EtaDetector",
     "The maximum allowed pseudorapidity for jets to be considered",
     &Cutter::theEtaDetector, 5.0, 0.0, 10.0, true, false, true);

  static Parameter<Cutter,double> interfacejetdefpt
    ("jetdefpt",
     "The minimum transversal momentum to define a jet",
     &Cutter::theJetDefPt, 15.0, 0.0, 14000, true, false, true );
}
