// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZAnalysis class.
//

#include "ZAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

ZAnalysis::ZAnalysis() : ZmassHadron_(82.,102.,100), ZmassDetector_(82.,102.,100)
{}

void ZAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  PGSInterface::analyze(event, ieve, loop, state);
  // detector level
  // require two opposite sign same flavour leptons
  unsigned int neminus(0),neplus(0),nmuminus(0),nmuplus(0);
  Lorentz5Momentum pz;
  for(unsigned int ix=0;ix<reconstructedObjects().size();++ix) {
    if(reconstructedObjects()[ix].type==Electron) {
      pz += reconstructedObjects()[ix].momentum;
      if(reconstructedObjects()[ix].PDGcode<0) ++neplus;
      else                                     ++neminus;
    }
    else if(reconstructedObjects()[ix].type==Muon) {
      pz += reconstructedObjects()[ix].momentum;
      if(reconstructedObjects()[ix].PDGcode<0) ++nmuplus;
      else                                     ++nmuminus;
    }
  }
  // require 2 leptons
  if(neminus+neplus+nmuminus+nmuplus!=2) return;
  // require opposite sign
  if(neminus!=neplus||nmuminus!=nmuplus) return;
  double mz = pz.m()/GeV;
  ZmassDetector_ += mz;
  // hadron level
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator stest =event->primaryCollision()->steps().end();
  StepVector::const_iterator send=sit;
  ++send;
  if(send==stest) --send;
  ++send;
  if(send==stest) --send;
  ++send;
  pz = LorentzMomentum();
  for(;sit!=send;++sit) {
    ParticleSet part=(**sit).all();
    ParticleSet::const_iterator iter = part.begin(), end = part.end();
    for( ;iter!=end;++iter) {
      if((**iter).children().size()!=2) continue;
      if((**iter).id()==ParticleID::Z0||(**iter).id()==ParticleID::gamma) {
	pz=(*iter)->momentum();
	double mz = pz.mass()/GeV;
	ZmassHadron_ += mz;
      }  
    }
  }
}

IBPtr ZAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr ZAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


NoPIOClassDescription<ZAnalysis> ZAnalysis::initZAnalysis;
// Definition of the static class description member.

void ZAnalysis::Init() {

  static ClassDocumentation<ZAnalysis> documentation
    ("There is no documentation for the ZAnalysis class");

}

void ZAnalysis::dofinish() {
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  string title;
  using namespace HistogramOptions;
  ZmassHadron_  .topdrawOutput(outfile,Frame,"BLACK","Z mass");
  ZmassDetector_.topdrawOutput(outfile,HistogramOptions::None,"RED");
  outfile.close();
  PGSInterface::dofinish();
}
