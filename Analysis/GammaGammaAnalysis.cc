// -*- C++ -*-
//
// GammaGammaAnalysis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGammaAnalysis class.
//

#include "GammaGammaAnalysis.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
GammaGammaAnalysis::GammaGammaAnalysis() :
  _ptharder(0.,250.,100), _ptsofter(0.,250.,100), _ptpair(0.,250.,100), 
  _Eharder(0.,3000.,100), _Esofter(0.,3000.,100), _Epair(0.,6000.,100), 
  _rapharder(-12.,12.,120), _rapsofter(-12.,12.,120), _rappair(-12.,12.,120), 
  _phiharder(-Constants::pi,Constants::pi,50), 
  _phisofter(-Constants::pi,Constants::pi,50), 
  _deltaphi(-Constants::pi,Constants::pi,100), 
  _mpair(0,1000,100)  {}

void GammaGammaAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  _ptharder.topdrawOutput(outfile,Frame|Ylog,"BLACK","pt harder");
  _ptsofter.topdrawOutput(outfile,Frame|Ylog,"BLACK","pt softer");
  _ptpair.topdrawOutput(outfile,Frame|Ylog,"BLACK","pt pair");
  _Eharder.topdrawOutput(outfile,Frame|Ylog,"BLACK","E harder");
  _Esofter.topdrawOutput(outfile,Frame|Ylog,"BLACK","E softer");
  _Epair.topdrawOutput(outfile,Frame|Ylog,"BLACK","E pair");
  _rapharder.topdrawOutput(outfile,Frame,"BLACK","y harder");
  _rapsofter.topdrawOutput(outfile,Frame,"BLACK","y softer");
  _rappair.topdrawOutput(outfile,Frame,"BLACK","y pair");
  _phiharder.topdrawOutput(outfile,Frame,"BLACK","phi harder");
  _phisofter.topdrawOutput(outfile,Frame,"BLACK","phi softer");
  _deltaphi.topdrawOutput(outfile,Frame,"BLACK","Delta phi");
  _mpair.topdrawOutput(outfile,Frame|Ylog,"BLACK","M pair");
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


void GammaGammaAnalysis::analyze(tEventPtr event, long, int, int) {
  //  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // find the Z
  Lorentz5Momentum p1, p2, ppair;  
  bool foundphotons = false;
  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));

  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if((**it).id()==ParticleID::gamma) {
      // find the two hardest photons in the event
      if( getMomentum(*it).perp() > p2.perp() ) {
	if (getMomentum(*it).perp() > p1.perp()) {
	  p2 = p1;
	  p1 = getMomentum(*it);
	} else {
	  p2 = getMomentum(*it);
	}
      }
    }
  }

  //  cerr << "E1 = " <<  p1.e()/GeV << ", E1 = " <<  p1.e()/GeV << "\n";
  if (p1.perp()/GeV > 0 && p2.perp()/GeV > 0) foundphotons = true;

  ppair = p1 + p2;
  if (foundphotons) {
    _ptharder += p1.perp()/GeV;
    _ptsofter += p2.perp()/GeV;
    _ptpair += ppair.perp()/GeV;
    _Eharder += p1.e()/GeV;
    _Esofter += p2.e()/GeV;
    _Epair += ppair.e()/GeV;
    _rapharder += p1.rapidity();
    _rapsofter += p2.rapidity();
    _rappair += ppair.rapidity();
    _phiharder += p1.phi();
    _phisofter += p2.phi();
    _deltaphi += (p2.vect()).deltaPhi(p1.vect());
    _mpair += ppair.m()/GeV;
  } else {
    cerr << "Analysis/GammaGammaAnalysis: Found no hard photon in event " 
	 << event->number()  << ".\n";
    generator()->log() << "Analysis/GammaGammaAnalysis: " 
		       << "Found no hard photon in event " 
		       << event->number()  << ".\n"
		       << *event;    
  }  
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<GammaGammaAnalysis,AnalysisHandler>
describeHerwigGammaGammaAnalysis("Herwig::GammaGammaAnalysis", "HwAnalysis.so");

void GammaGammaAnalysis::Init() {

  static ClassDocumentation<GammaGammaAnalysis> documentation
    ("There is no documentation for the GammaGammaAnalysis class");

}

