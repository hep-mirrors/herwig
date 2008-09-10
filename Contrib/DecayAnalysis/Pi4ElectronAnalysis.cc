// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Pi4ElectronAnalysis class.
//

#include "Pi4ElectronAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include <ThePEG/EventRecord/Event.h>
#include <ThePEG/PDT/EnumParticles.h>
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

void Pi4ElectronAnalysis::analyze(tEventPtr event, long , int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  // find all the decaying neutral pions
  tPVector particles;
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
    ThePEG::ParticleSet::iterator iter=part.begin();
    ThePEG::ParticleSet::iterator end=part.end();
    for( ;iter!=end;++iter) {
      if((**iter).id()==ParticleID::pi0) {
	particles.push_back(*iter);
      }
    }
  }
  // analyse them
  analyze(particles);
}

void Pi4ElectronAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void Pi4ElectronAnalysis::analyze(tPPtr part) {
  // check the number of children of the particle
  if(part->children().size()!=4)  return;
  Lorentz5Momentum pe[2],pp[2];
  unsigned int ne(0),np(0);
  int id;
  // find the particles
  unsigned int ix,iy;
  for(ix=0;ix<part->children().size();++ix) {
    id=part->children()[ix]->id();
    if(id==ParticleID::eplus){pe[ne]=part->children()[ix]->momentum();++ne;}
    else if(id==ParticleID::eminus){pp[np]=part->children()[ix]->momentum();++np;}
  }
  if(ne!=2||np!=2){return;}
  // find the invariant masses
  Lorentz5Momentum ptemp;
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      ptemp=pe[ix]+pp[iy];
      ptemp.rescaleMass();
      for(unsigned int iz=0;iz<_mffbar.size();++iz) {
	*_mffbar[iz]+=ptemp.mass()/MeV;
      }
    }
  }
}

NoPIOClassDescription<Pi4ElectronAnalysis> Pi4ElectronAnalysis::initPi4ElectronAnalysis;
// Definition of the static class description member.

void Pi4ElectronAnalysis::Init() {

  static ClassDocumentation<Pi4ElectronAnalysis> documentation
    ("There is no documentation for the Pi4ElectronAnalysis class");

}

void Pi4ElectronAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  for(unsigned int ix=0;ix<_mffbar.size();++ix) {
    _mffbar[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "Mass of the e2+3e2-3 pair in P203Re2+3e2-3e2+3e2-3",
			       "             X X X X         GX XW X X X X X X X X",
			       "1/SdS/dm0e2+3e2-31/GeV2-13",
			       "  G G   X X X X XX    X  X",
			       "m0e2+3e2-31",
			       " X X X X XX");
  }
}

void Pi4ElectronAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  _mffbar.push_back(new_ptr(Histogram(0.0,140.0,1000)));
  _mffbar.push_back(new_ptr(Histogram(0.0, 20.0,100 )));
}
