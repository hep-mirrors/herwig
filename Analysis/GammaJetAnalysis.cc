// -*- C++ -*-
//
// GammaJetAnalysis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaJetAnalysis class.
//

#include "GammaJetAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

GammaJetAnalysis::GammaJetAnalysis() :
  _ptg(0.,250.,100), _ptgZoom(35.,65.,100), 
  _Eg(0,3000,100), _rapg(-10.,10.,100), 
  _phig(-Constants::pi,Constants::pi,100) {}

void GammaJetAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  _Eg.topdrawOutput(outfile,Frame,"BLACK","Energy of Gamma");
  _Eg.topdrawOutput(outfile,Frame|Ylog,"BLACK","Energy of Gamma");
  _ptg.topdrawOutput(outfile,Frame,"BLACK","pt of Gamma");
  _ptg.topdrawOutput(outfile,Frame|Ylog,"BLACK","pt of Gamma");
  _ptgZoom.topdrawOutput(outfile,Frame,"BLACK","35<pt/GeV<65 of Gamma");
  _ptgZoom.topdrawOutput(outfile,Frame|Ylog,"BLACK","35<pt/GeV<65 of Gamma");
  _rapg.topdrawOutput(outfile,Frame,"BLACK","Rapidity of Gamma");
  _rapg.topdrawOutput(outfile,Frame|Ylog,"BLACK","Rapidity of Gamma");
  _phig.topdrawOutput(outfile,Frame,"BLACK","Azimuth of Gamma");
  _phig.topdrawOutput(outfile,Frame|Ylog,"BLACK","Azimuth of Gamma");
  outfile.close();
}

namespace {
  inline Lorentz5Momentum getMomentum(tcPPtr particle) {
    return particle->momentum();
    //Lorentz5Momentum tmp = particle->children()[0]->next()->momentum();
    //tmp += particle->children()[1]->next()->momentum();
    //tmp.rescaleMass();
    //return tmp;

  }
}


void GammaJetAnalysis::analyze(tEventPtr event, long, int, int) {
  //  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // find the Z
  Lorentz5Momentum pg;  
  bool foundphoton = false;
  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));

  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if((**it).id()==ParticleID::gamma) {
      // only book the hardest photon in the event
      if( (**it).momentum().perp() > pg.perp() ) {
	foundphoton = true;
	pg=getMomentum(*it);
      }
    }
  }

  if (foundphoton) {
    Energy pt = pg.perp();
    (_ptg)+=(pt)/GeV;
    (_Eg)+=pg.e()/GeV;
    (_ptgZoom)+=(pt)/GeV;
    double rap = 0.5*log((pg.e()+pg.z())/(pg.e()-pg.z()));
    (_rapg)+=(rap);
    (_phig)+=pg.phi();
  } else {
    cerr << "Analysis/GammaJetAnalysis: Found no hard photon in event " 
	 << event->number()  << ".\n";
    generator()->log() << "Analysis/GammaJetAnalysis: " 
		       << "Found no hard photon in event " 
		       << event->number()  << ".\n"
		       << *event;    
  }  
}

NoPIOClassDescription<GammaJetAnalysis> GammaJetAnalysis::initGammaJetAnalysis;
// Definition of the static class description member.

void GammaJetAnalysis::Init() {

  static ClassDocumentation<GammaJetAnalysis> documentation
    ("There is no documentation for the GammaJetAnalysis class");

}

