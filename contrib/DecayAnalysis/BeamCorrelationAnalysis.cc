// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BeamCorrelationAnalysis class.
//

#include "BeamCorrelationAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;

void BeamCorrelationAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  PPair beam=event->primaryCollision()->primarySubProcess()->incoming();
  ParticleVector incoming=event->primaryCollision()->primarySubProcess()->outgoing();
  if(incoming.size()!=1) return;
  ParticleVector outgoing=incoming[0]->children();
  if(outgoing.size()!=2) return;
  int id0=incoming[0]->id();
  int id1=abs(outgoing[0]->id());
  int iloc(-1);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(id0==_incoming[ix]&&id1==_outgoing[ix]) {
      iloc=ix;
      break;
    }
  }
  if(iloc<0) {
    iloc=_incoming.size();
    _incoming.push_back(id0);
    _outgoing.push_back(id1);
    _angle.push_back(new_ptr(Histogram(Histogram(0.,Constants::pi,200))));
  }
  tPPtr in = beam.first->id()>0 ? beam.first : beam.second;
  tPPtr out;
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    if(outgoing[ix]->id()>0) out=outgoing[ix];
  }
  if(!in||!out) return;
  _angle[iloc]->addWeighted(in->momentum().angle(out->momentum()),
			    event->weight());
}

NoPIOClassDescription<BeamCorrelationAnalysis> 
BeamCorrelationAnalysis::initBeamCorrelationAnalysis;
// Definition of the static class description member.

void BeamCorrelationAnalysis::Init() {

  static ClassDocumentation<BeamCorrelationAnalysis> documentation
    ("There is no documentation for the BeamCorrelationAnalysis class");

}

inline void BeamCorrelationAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + 
    string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    string title= "Angle for " + getParticleData(_incoming[ix])->PDGName() +
      " -> " + getParticleData(_outgoing[ix])->PDGName() +
      getParticleData(-_outgoing[ix])->PDGName();
    using namespace HistogramOptions;
    _angle[ix]->topdrawOutput(output,Frame|Errorbars,
			      "RED",title,"",
			      "1/SdS/dQ",
			      "  G G  G",
			      "Q",
			      "G");
  }
}
