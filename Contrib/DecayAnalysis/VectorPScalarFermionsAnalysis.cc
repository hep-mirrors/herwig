// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorPScalarFermionsAnalysis class.
//

#include "VectorPScalarFermionsAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

void VectorPScalarFermionsAnalysis::analyze(tEventPtr event, long, int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  // find all vector particles 
  tPVector particles;
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
    ThePEG::ParticleSet::iterator iter=part.begin();
    ThePEG::ParticleSet::iterator end=part.end();
    for( ;iter!=end;++iter) {
      if((**iter).dataPtr()->iSpin()==PDT::Spin1&&(**iter).id()!=ParticleID::gamma&&
	 (**iter).children().size()==3) {
	if((**iter).children()[0]->dataPtr()->iSpin()==PDT::Spin0)
	  particles.push_back(*iter);
      }
    }
  }
  // analyse them
  analyze(particles);
}

void VectorPScalarFermionsAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void VectorPScalarFermionsAnalysis::analyze(tPPtr part) {
  long id[3]={part->children()[0]->id(),part->children()[1]->id(),
	     part->children()[2]->id()};
  // id's of the fermions
  if(abs(id[1])!=abs(id[2])){return;}
  if(abs(id[1])!=11&&abs(id[1])!=13){return;}
  // check if we already have this decay
  unsigned int ix=0; bool found(false); 
  while(!found&&ix<_incoming.size()) {
    if(_incoming[ix]==part->id()&&_outgoingP[ix]==id[0]&&_outgoingf[ix]==abs(id[1])) {
      found=true;
    }
    else {
      ++ix;
    }
  }  // create a new graph if needed
  if(!found) {
    ix=_incoming.size();
    _incoming.push_back(part->id());
    _outgoingP.push_back(id[0]);
    _outgoingf.push_back(abs(id[1]));
    _mffa.push_back(new_ptr(Histogram(0.0,
				     (part->nominalMass()+part->dataPtr()->widthUpCut())/MeV,
				     200)));
    _mffb.push_back(new_ptr(Histogram(0.0,20,200)));
    _mPf.push_back(new_ptr(Histogram(0.0,
				     (part->nominalMass()+part->dataPtr()->widthUpCut())/MeV,
				     200)));
    _mPfbar.push_back(new_ptr(Histogram(0.0,
				     (part->nominalMass()+part->dataPtr()->widthUpCut())/MeV,
				     200)));
  }
  // add the results to the histogram
  Lorentz5Momentum ptemp;
  ptemp=part->children()[1]->momentum()+part->children()[2]->momentum();
  ptemp.rescaleMass();
  *_mffa[ix]+=ptemp.mass()/MeV;
  *_mffb[ix]+=ptemp.mass()/MeV;
  ptemp=part->children()[0]->momentum()+part->children()[1]->momentum();
  ptemp.rescaleMass();
  *_mPf[ix]+=ptemp.mass()/MeV;
  ptemp=part->children()[0]->momentum()+part->children()[2]->momentum();
  ptemp.rescaleMass();
  *_mPfbar[ix]+=ptemp.mass()/MeV;
}

NoPIOClassDescription<VectorPScalarFermionsAnalysis> VectorPScalarFermionsAnalysis::initVectorPScalarFermionsAnalysis;
// Definition of the static class description member.

void VectorPScalarFermionsAnalysis::Init() {

  static ClassDocumentation<VectorPScalarFermionsAnalysis> documentation
    ("There is no documentation for the VectorPScalarFermionsAnalysis class");

}


void VectorPScalarFermionsAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  string title,temp;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    title= getParticleData(_incoming[ix])->PDGName() +
      " -> " + getParticleData(_outgoingP[ix])->PDGName()+  " " +
      getParticleData(_outgoingf[ix])->PDGName() + " " +
      getParticleData(-_outgoingf[ix])->PDGName();
    temp = "Mass for f fbar in " +title;
    using namespace HistogramOptions;
    _mffa[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",temp,"",
			     "1/SdS/dm0l2+3l2-31/GeV2-13",
			     "  G G   X X X X XX    X  X",
			     "m0l2+3l2-31",
			     " X X X X XX");
    _mffb[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",temp,"",
			     "1/SdS/dm0l2+3l2-31/GeV2-13",
			     "  G G   X X X X XX    X  X",
			     "m0l2+3l2-31",
			     " X X X X XX");
    temp = "Mass for vector fermion mass in " +title;
    _mPf[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",temp,"",
			    "1/SdS/dm0Pl2-31/GeV2-13",
			    "  G G   X  X XX    X  X",
			    "m0Pl2-31",
			    " X  X XX");
    temp = "Mass for vector fbar mass    in " +title;
    _mPfbar[ix]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",temp,"",
			       "1/SdS/dm0Pl2+31/GeV2-13",
			       "  G G   X  X XX    X  X",
			       "m0Pl2+31",
			       " X  X XX");
  }
}
