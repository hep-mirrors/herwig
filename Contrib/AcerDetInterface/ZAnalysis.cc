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
  AcerDet::analyze(event, ieve, loop, state);
  // detector level
  // require 2 leptons
  if(numberOfLeptons()!=2) return;
  // require opposite sign
  if(leptonID()[0]!=-leptonID()[1]) return;
  Lorentz5Momentum pz = leptonMomentum()[0]+leptonMomentum()[1];
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
  ZmassDetector_.topdrawOutput(outfile,None,"RED");
  outfile.close();
  AcerDet::dofinish();
}
