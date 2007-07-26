// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the a1DecayAnalysis class.
//

#include "a1DecayAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;

void a1DecayAnalysis::analyze(tEventPtr event, long, int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  // find all omega and phi particles 
  tPVector particles;
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
    ThePEG::ParticleSet::iterator iter=part.begin();
    ThePEG::ParticleSet::iterator end=part.end();
    for( ;iter!=end;++iter) {
      int id=abs((**iter).id());
      if(id==ParticleID::a_10||id==ParticleID::a_1plus) 
	particles.push_back(*iter);
    }
  }
  // analyse them
  analyze(particles);
}

void a1DecayAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void a1DecayAnalysis::analyze(tPPtr part) {
  // find the pions
  ParticleVector pions;
  findPions(part,pions);
  if(pions.size()!=3) return;
}

NoPIOClassDescription<a1DecayAnalysis> a1DecayAnalysis::inita1DecayAnalysis;
// Definition of the static class description member.

void a1DecayAnalysis::Init() {

  static ClassDocumentation<a1DecayAnalysis> documentation
    ("There is no documentation for the a1DecayAnalysis class");

}

void a1DecayAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = CurrentGenerator::current().filename() + 
    string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  _hist0->topdrawOutput(output,true,true,false,false,"RED",
			"P203P203 mass in A011203RP203P203P203",
			"GX XGX X          X XX XWGX XGX XGX X",
			"1/GdG/dm0P203P2031/MeV2-13",
			"  F F   XGX XGX XX    X  X",
			"m0P203P2031/MeV",
			" XGX XGX XX    ");
  _hist1A->topdrawOutput(output,true,true,false,false,"RED",
			 "P203P203 mass in A0112+3RP203P203P2+3",
			 "GX XGX X          X XX XWGX XGX XGX X",
			 "1/GdG/dm0P203P2031/MeV2-13",
			 "  F F   XGX XGX XX    X  X",
			 "m0P203P2031/MeV",
			 " XGX XGX XX    ");
  _hist1B->topdrawOutput(output,true,true,false,false,"RED",
			 "P203P203 mass in A0112+3RP203P203P2+3",
			 "GX XGX X          X XX XWGX XGX XGX X",
			 "1/GdG/dm0P203P2+31/MeV2-13",
			 "  F F   XGX XGX XX    X  X",
			 "m0P203P2+31/MeV",
			 " XGX XGX XX    ");
  _hist2A->topdrawOutput(output,true,true,false,false,"RED",
			"P203P203 mass in A011203RP2+3P2-3P203",
			"GX XGX X          X XX XWGX XGX XGX X",
			"1/GdG/dm0P2+3P2-31/MeV2-13",
			"  F F   XGX XGX XX    X  X",
			"m0P2+3P2-31/MeV",
			" XGX XGX XX    ");
  _hist2B->topdrawOutput(output,true,true,false,false,"RED",
			"P203P203 mass in A011203RP2+3P2-3P203",
			"GX XGX X          X XX XWGX XGX XGX X",
			"1/GdG/dm0P2+3P2031/MeV2-13",
			"  F F   XGX XGX XX    X  X",
			"m0P2+3P2031/MeV",
			" XGX XGX XX    ");
  _hist2C->topdrawOutput(output,true,true,false,false,"RED",
			"P203P203 mass in A011203RP2+3P2-3P203",
			"GX XGX X          X XX XWGX XGX XGX X",
			"1/GdG/dm0P2-3P2031/MeV2-13",
			"  F F   XGX XGX XX    X  X",
			"m0P2-3P2031/MeV",
			" XGX XGX XX    ");
  _hist3A->topdrawOutput(output,true,true,false,false,"RED",
			 "P203P203 mass in A0112+3RP2+3P2+3P2-3",
			 "GX XGX X          X XX XWGX XGX XGX X",
			 "1/GdG/dm0P2+3P2+31/MeV2-13",
			 "  F F   XGX XGX XX    X  X",
			 "m0P2+3P2+31/MeV",
			 " XGX XGX XX    ");
  _hist3A->topdrawOutput(output,true,true,false,false,"RED",
			 "P203P203 mass in A0112+3RP2+3P2+3P2-3",
			 "GX XGX X          X XX XWGX XGX XGX X",
			 "1/GdG/dm0P2+3P2-31/MeV2-13",
			 "  F F   XGX XGX XX    X  X",
			 "m0P2+3P2-31/MeV",
			 " XGX XGX XX    ");
}

void a1DecayAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  _hist0  = new_ptr(Histogram(0.,1500.,200));
  _hist1A = new_ptr(Histogram(0.,1500.,200));
  _hist1B = new_ptr(Histogram(0.,1500.,200));
  _hist2A = new_ptr(Histogram(0.,1500.,200));
  _hist2B = new_ptr(Histogram(0.,1500.,200));
  _hist2C = new_ptr(Histogram(0.,1500.,200));
  _hist3A = new_ptr(Histogram(0.,1500.,200));
  _hist3B = new_ptr(Histogram(0.,1500.,200));
}

void a1DecayAnalysis::findPions(tPPtr part,ParticleVector & pions) {
  if(abs(part->id())==ParticleID::piplus||part->id()==ParticleID::pi0) {
    pions.push_back(part);
    return;
  }
  else if(!part->children().empty()) {
    for(unsigned int ix=0;ix<part->children().size();++ix) {
      findPions(part->children()[ix],pions);
    }
  }
}
