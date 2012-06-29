// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoSGammaAnalysis class.
//

#include "BtoSGammaAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

void BtoSGammaAnalysis::analyze(tEventPtr event, long, int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  // find all three body B decays
  tPVector particles;
  for(unsigned int ix=1, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part(event->primaryCollision()->step(ix)->all());
    ThePEG::ParticleSet::iterator iter(part.begin()),end(part.end());
    for( ;iter!=end;++iter) {
      if((abs((**iter).id())==ParticleID::Bplus||abs((**iter).id())==ParticleID::B0))
	  particles.push_back(*iter);
    }
  }
  // analyse them
  analyze(particles);
}

void BtoSGammaAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void BtoSGammaAnalysis::analyze(tPPtr part) {
  Lorentz5Momentum phad;
  Lorentz5Momentum pgamma;
  unsigned int ngamma(0);
  for(unsigned int ix=0;ix<part->children().size();++ix) {
    if(part->children()[ix]->id()==ParticleID::gamma) {
      ++ngamma;
      pgamma+=part->children()[ix]->momentum();
    }
    else {
      phad+=part->children()[ix]->momentum();
    }
  }
  if(ngamma!=1) return;
  int id(part->id());
  unsigned int imode;
  if(id==ParticleID::B0)          imode=0;
  else if(id==ParticleID::Bbar0)  imode=1;
  else if(id==ParticleID::Bplus)  imode=2;
  else if(id==ParticleID::Bminus) imode=3;
  else                            return;
  // calculate the hadronic mass
  phad.rescaleMass();
  *_hadmass[imode] += phad.mass()/MeV;
  pgamma.boost(-part->momentum().boostVector());
  *_spectrum[imode]+=pgamma.e()/MeV;
}

NoPIOClassDescription<BtoSGammaAnalysis> BtoSGammaAnalysis::initBtoSGammaAnalysis;
// Definition of the static class description member.

void BtoSGammaAnalysis::Init() {

  static ClassDocumentation<BtoSGammaAnalysis> documentation
    ("There is no documentation for the BtoSGammaAnalysis class");

}

void BtoSGammaAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() 
    + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  // output the histograms
  string title,temp[2],tcase;
  for(unsigned int ix=0;ix<4;++ix) {
    if(ix==0) {
      title="B203";
      tcase=" X X";
    }
    else if(ix==1) {
      title="B0O203";
      tcase=" UDX X";
    }
    else if(ix==2) {
      title="B2+3";
      tcase=" X X";
    }
    else if(ix==3) {
      title="B2-3";
      tcase=" X X";
    }
    temp[0] = "Hadronic Mass for " + title;
    temp[1] = "                  " + tcase;
    using namespace HistogramOptions;
    _hadmass[ix]->topdrawOutput(output,Frame|Errorbars,
				"RED",
				temp[0],temp[1],
				"1/SdS/dm/GeV2-13",
				"  G G       X  X",
				"m/GeV",
				"     ");
    temp[0]= "Photon spectrum for " + title;
    temp[1]= "                    " + tcase;
    _spectrum[ix]->topdrawOutput(output,Frame|Errorbars,
				 "RED",
				 temp[0],temp[1],
				 "1/SdS/dE0G1/GeV2-13",
				 "  G G   XGX    X  X",
				 "E0G1/GeV",
				 " XGX    ");
  }
}

void BtoSGammaAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  for(unsigned int ix=0;ix<4;++ix) {
    _hadmass.push_back(new_ptr(Histogram(850.,4500.,100)));
    _spectrum.push_back(new_ptr(Histogram(0.,2800.,100)));
  }
}
