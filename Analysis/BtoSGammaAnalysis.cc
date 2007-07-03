// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoSGammaAnalysis class.
//

#include "BtoSGammaAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"


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
    for( ;iter!=end;++iter)
      {if((abs((**iter).id())==ParticleID::Bplus||abs((**iter).id())==ParticleID::B0)
	  &&(**iter).children().size()==3){particles.push_back(*iter);}}
  }
  // analyse them
  analyze(particles);
}

void BtoSGammaAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void BtoSGammaAnalysis::analyze(tPPtr part) {
  // check we have the right decay
  // first daugter is a strange
  if(abs(part->children()[0]->id())!=ParticleID::s) return;
  // second is the photon
  if(abs(part->children()[1]->id())!=ParticleID::gamma) return;
  // last is the fdinal quarl
  if(abs(part->children()[2]->id())>ParticleID::u) return;
  // work out the type of B
  int id(part->id());
  unsigned int imode;
  if(id==ParticleID::B0)          imode=0;
  else if(id==ParticleID::Bbar0)  imode=1;
  else if(id==ParticleID::Bplus)  imode=2;
  else if(id==ParticleID::Bminus) imode=3;
  else                            return;
  // calculate the hadronic mass
  Lorentz5Momentum ptemp;
  ptemp = part->children()[0]->momentum()+part->children()[1]->momentum();
  ptemp.rescaleMass();
  *_hadmass[imode] +=ptemp.mass()/MeV;
  ptemp= part->children()[1]->momentum();
  ptemp.boost(-part->momentum().boostVector());
  *_spectrum[imode]+=ptemp.e()/MeV;
}

NoPIOClassDescription<BtoSGammaAnalysis> BtoSGammaAnalysis::initBtoSGammaAnalysis;
// Definition of the static class description member.

void BtoSGammaAnalysis::Init() {

  static ClassDocumentation<BtoSGammaAnalysis> documentation
    ("There is no documentation for the BtoSGammaAnalysis class");

}

void BtoSGammaAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  ofstream output("BtoSGamma.top");
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
      _hadmass[ix]->topdrawOutput(output,true,true,false,true,
				  "RED",
				  temp[0],temp[1],
				  "1/SdS/dm/GeV2-13",
				  "  G G       X  X",
				  "m/GeV",
				  "     ");
      temp[0]= "Photon spectrum for " + title;
      temp[1]= "                    " + tcase;
      _spectrum[ix]->topdrawOutput(output,true,true,false,true,
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
    _hadmass.push_back(new_ptr(Histogram(0.,5300.,100)));
    _spectrum.push_back(new_ptr(Histogram(0.,5300.,100)));
  }
}
