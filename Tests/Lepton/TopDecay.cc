// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TopDecay class.
//

#include "TopDecay.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Herwig;

void TopDecay::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tPVector particles;
  event->selectFinalState(back_inserter(particles));
  vector<fastjet::PseudoJet> fastjet_particles;
  for (unsigned int j=0; j<particles.size(); j++) {
    tPPtr parent=particles[j];
    do {
      if(parent->parents().empty()) parent = tPPtr();
      else parent = parent->parents()[0];
    }
    while (parent&&abs(parent->id())!=ParticleID::Wplus);
    if(parent&&abs(parent->id())==ParticleID::Wplus) continue;
    fastjet::PseudoJet p(particles[j]->momentum().x()/GeV, 
			 particles[j]->momentum().y()/GeV, 
			 particles[j]->momentum().z()/GeV, 
			 particles[j]->momentum().e()/GeV);
    p.set_user_index(j);
    fastjet_particles.push_back(p);
  }
  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::JetDefinition jet_def(fastjet::ee_kt_algorithm, 
				 recomb_scheme, strategy);
  fastjet::ClusterSequence cs(fastjet_particles, jet_def);
  
  // ynm distributions
  y3_ += log10(cs.exclusive_ymerge(3));
}


IBPtr TopDecay::clone() const {
  return new_ptr(*this);
}

IBPtr TopDecay::fullclone() const {
  return new_ptr(*this);
}
void TopDecay::persistentOutput(PersistentOStream & os) const {
}

void TopDecay::persistentInput(PersistentIStream & is, int) {
}

ClassDescription<TopDecay> TopDecay::initTopDecay;
// Definition of the static class description member.

void TopDecay::Init() {

  static ClassDocumentation<TopDecay> documentation
    ("There is no documentation for the TopDecay class");

}

void TopDecay::dofinish() {
  AnalysisHandler::dofinish();
  using namespace HistogramOptions;
  string fname = generator()->filename() + string("-") 
    + name() + string(".top");
  ofstream output(fname.c_str());
  y3_.topdrawOutput(output,Frame,
		    "RED",
		    "y0231 ",
		    " X  X ",
		    "1/NdN/dlog(y0231)",
		    "            X  X ",
		    "log(y0231)",
		    "     X  X ");
}
