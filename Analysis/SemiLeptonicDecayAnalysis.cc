// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SemiLeptonicDecayAnalysis class.
//

#include "SemiLeptonicDecayAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;

void SemiLeptonicDecayAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {  
  // Rotate to CMS, extract final state particles and call analyze(particles).
  AnalysisHandler::analyze(event, ieve, loop, state);
  tPVector particles;
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
    ThePEG::ParticleSet::iterator iter=part.begin();
    ThePEG::ParticleSet::iterator end=part.end();
    for( ;iter!=end;++iter) {
      if((**iter).dataPtr()->iSpin()==PDT::Spin0&&
	 (**iter).children().size()>=2&&(**iter).children().size()<=3) {
	particles.push_back(*iter);
      }
    }
  }
  // analyse them
  analyze(particles);
}

void SemiLeptonicDecayAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void SemiLeptonicDecayAnalysis::analyze(tPPtr part) {
  // check the number of children of the particle
  // check this is a decay
  tPPtr lep[2];
  if(part->children().size()==2&&
     abs(part->children()[1]->id())==24&&part->children()[1]->children().size()==2) {
      lep[0]=part->children()[1]->children()[0];
      lep[1]=part->children()[1]->children()[1];
    }
  else if(part->children().size()==3&&
	  abs(part->children()[1]->id())>=11&&abs(part->children()[1]->id())<=16&&
	  abs(part->children()[2]->id())>=11&&abs(part->children()[2]->id())<=16) {
      lep[0]=part->children()[1];
      lep[1]=part->children()[2];
    }
  else {
    return;
  }
  // find the ids of the lepton and the neutrino
  int id1(abs(lep[0]->id())),id2(abs(lep[1]->id()));
  unsigned int ilep(0),inu(0),loce(1),locn(1);
  if(id1%2==0&&id1>11&&id1<17){inu=(id1-10)/2;locn=0;}
  else if(id2%2==0&&id2>11&&id2<17){inu=(id2-10)/2;locn=1;}
  if(id1%2==1&&id1>10&&id1<16){ilep=(id1-9)/2;loce=0;}
  else if(id2%2==1&&id2>10&&id2<16){ilep=(id2-9)/2;loce=1;}
  if(ilep==0||inu==0){return;}
  unsigned int ix=0; bool found(false);
  while(!found&&ix<_incoming.size()) {
    if(_incoming[ix]==part->id()&&_outgoing[ix]==part->children()[0]->id()&&
       ilep==_outgoingL[ix]){found=true;}
    else{++ix;}
  }
  if(!found) {
    ix=_incoming.size();
    _incoming.push_back(part->id());
    _outgoing.push_back(part->children()[0]->id());
    _outgoingL.push_back(ilep);
    _energy.push_back(new_ptr(Histogram(0.0,
					(part->mass()-part->children()[0]->mass())/MeV,
					200)));
    _scale.push_back(new_ptr(Histogram(0.0,
				       (part->mass()-part->children()[0]->mass())/MeV,
				       200)));
  }
  // add the results to the histogram
  Lorentz5Momentum ptemp;
  ptemp = lep[0]->momentum()+lep[1]->momentum();
  ptemp.rescaleMass();
  *_scale[ix] += ptemp.mass()/MeV;
  ptemp = part->children()[0]->momentum()+lep[locn]->momentum();
  ptemp.rescaleMass();
  Energy ee = 1./2./part->mass()*
    ( part->mass()*part->mass()-ptemp.mass()*ptemp.mass()
      +lep[loce]->mass()*lep[loce]->mass());
  *_energy[ix] += ee/MeV;
}

NoPIOClassDescription<SemiLeptonicDecayAnalysis> SemiLeptonicDecayAnalysis::initSemiLeptonicDecayAnalysis;
// Definition of the static class description member.

void SemiLeptonicDecayAnalysis::Init() {

  static ClassDocumentation<SemiLeptonicDecayAnalysis> documentation
    ("There is no documentation for the SemiLeptonicDecayAnalysis class");

}

void SemiLeptonicDecayAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  string title,temp;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    title= getParticleData(_incoming[ix])->PDGName() +
      " -> " + getParticleData(_outgoing[ix])->PDGName()+  " " +
      getParticleData(9+2*_outgoingL[ix])->PDGName() + " " +
      getParticleData(10+2*_outgoingL[ix])->PDGName();
    temp = "Mass for l nu in " +title;
    using namespace HistogramOptions;
    _scale[ix]->topdrawOutput(output,Frame|Errorbars,
			     "RED",temp,"",
			     "1/GdG/dm0lN1/MeV2-13",
			     "  F F   X GX    X  X",
			     "m0lN1/MeV",
			     " X GX    ");
    temp = "Lepton energy for in " +title;
    _energy[ix]->topdrawOutput(output,Frame|Errorbars,
			     "RED",temp,"",
			     "1/GdG/dE0l1/MeV2-13",
			     "  F F   X X    X  X",
			     "E0l1/MeV",
			     " X X    ");
  }
}

