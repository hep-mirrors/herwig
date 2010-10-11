// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LPairAnalysis class.
//

#include "LPairAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

namespace {
  inline Lorentz5Momentum getMomentum(tcPPtr particle) {
    return particle->momentum();
  }

  inline bool isLeptonPlus(tcPPtr p) {
    if ( p->id() == ParticleID::eplus ||
	 p->id() == ParticleID::muplus ||
	 p->id() == ParticleID::tauplus ) {
      return true;
    }else {
      return false;
    }
  }

  inline bool isLeptonMinus(tcPPtr p) {
    if ( p->id() == ParticleID::eminus ||
	 p->id() == ParticleID::muminus ||
	 p->id() == ParticleID::tauminus ) {
      return true;
    }else {
      return false;
    }
  }
  inline bool isFromTop(tcPPtr p) {
    while (p->parents()[0] && p->parents().size() == 1) {
      p = p->parents()[0];
      if (abs(p->id()) == ParticleID::t) {
	return true;
      }
    }
    return false; 
  } 
}

LPairAnalysis::LPairAnalysis() :
  _ptp(0.,250.,100), _ptm(0.,250.,100), _ptpair(0.,250.,100), 
  _etp(0.,250.,100), _etm(0.,250.,100), _etpair(0.,250.,100), 
  _ep(0.,1000.,100), _em(1.,3000.,100), _epair(0.,1500.,100), 
  _rapp(-5.,5.,100), _rapm(-5.,5.,100), _rappair(-5.,5.,100), 
  _phip(-Constants::pi,Constants::pi,50), 
  _phim(-Constants::pi,Constants::pi,50), 
  _deltaphi(-Constants::pi,Constants::pi,100), _mpair(0,750,100), 
  _etsum(0.,400.,100), _ptsum(0.,400.,100)
{}

void LPairAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  _ptp.topdrawOutput(outfile,Frame|Ylog,"RED","pt lp, lm");
  _ptm.topdrawOutput(outfile,Ylog,"BLUE","");
  _ptpair.topdrawOutput(outfile,Frame|Ylog,"BLACK","pt pair");
  _etp.topdrawOutput(outfile,Frame|Ylog,"RED","Et lp, lm");
  _etm.topdrawOutput(outfile,Ylog,"BLUE","");
  _etpair.topdrawOutput(outfile,Frame|Ylog,"BLACK","Et pair");
  _ep.topdrawOutput(outfile,Frame|Ylog,"RED","E lp, lm");
  _em.topdrawOutput(outfile,Ylog,"BLUE","");
  _epair.topdrawOutput(outfile,Frame|Ylog,"BLACK","E pair");
  _rapp.topdrawOutput(outfile,Frame,"RED","y lp, lm");
  _rapm.topdrawOutput(outfile,None,"BLUE","");
  _rappair.topdrawOutput(outfile,Frame,"BLACK","y pair");
  _phip.topdrawOutput(outfile,Frame,"RED","phi lp, lm");
  _phim.topdrawOutput(outfile,None,"BLUE","");
  _deltaphi.topdrawOutput(outfile,Frame,"BLACK","Delta phi");
  _mpair.topdrawOutput(outfile,Frame|Ylog,"BLACK","M pair");
  _etsum.topdrawOutput(outfile,Frame|Ylog,"BLACK","pt sum");
  _ptsum.topdrawOutput(outfile,Frame|Ylog,"BLACK","Et sum");
}

void LPairAnalysis::analyze(tEventPtr event, long, int, int) {
  Lorentz5Momentum ppair, plp, plm;  
  bool foundlp = false;
  bool foundlm = false;
  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));

  // find highest pt lepton+ and lepton- in the event resp.
  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if( isLeptonPlus(*it) ) {
      //     if ( getMomentum(*it).perp() > plp.perp() ) {
      if (isFromTop(*it)) {
	plp = getMomentum(*it);
	foundlp = true;
      }
    } else if( isLeptonMinus(*it) ) {
      //      if ( getMomentum(*it).perp() > plm.perp() ) {
      if (isFromTop(*it)) {
	plm = getMomentum(*it);
	foundlm = true;
      }
    }
  }
  
  if (foundlp && foundlm) {
    ppair = plp + plm;
    _ptp += plp.perp()/GeV;
    _ptm += plm.perp()/GeV;
    _ptpair += ppair.perp()/GeV;
    _etp += plp.et()/GeV;
    _etm += plm.et()/GeV;
    _etpair += ppair.et()/GeV;
    _ep += plp.e()/GeV;
    _em += plm.e()/GeV;
    _epair += ppair.e()/GeV;
    _rapp += plp.rapidity();
    _rapm += plm.rapidity();
    _rappair += ppair.rapidity();
    _phip += plp.phi();
    _phim += plm.phi();
    _deltaphi += (plp.vect()).deltaPhi(plm.vect());
    _mpair += ppair.m()/GeV;
    _etsum += (plp.et() + plm.et())/GeV;
    _ptsum += (plp.perp() + plm.perp())/GeV;
  } else {
    cerr << "Analysis/LPairAnalysis: did not find suitable lepton"
	 << " pair in event " << event->number()  << " ("
	 << (foundlp ? "+" : "0") 
	 << (foundlm ? "-" : "0") 
	 << ").\n";
    generator()->log() << "Analysis/LPairAnalysis: " 
		       << "Found no suitable lepton pair in event " 
		       << event->number()  << ".\n"
		       << *event;    
  }  
}

NoPIOClassDescription<LPairAnalysis> LPairAnalysis::initLPairAnalysis;
// Definition of the static class description member.

void LPairAnalysis::Init() {

  static ClassDocumentation<LPairAnalysis> documentation
    ("There is no documentation for the LPairAnalysis class");

}

