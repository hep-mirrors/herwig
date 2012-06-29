// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OniumToOniumPiPiAnalysis class.
//

#include "OniumToOniumPiPiAnalysis.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

void OniumToOniumPiPiAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  // Rotate to CMS, extract final state particles and call analyze(particles).
  AnalysisHandler::analyze(event, ieve, loop, state);
  tPVector particles;
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
    ThePEG::ParticleSet::iterator iter=part.begin();
    ThePEG::ParticleSet::iterator end=part.end();
    for( ;iter!=end;++iter) {
      if(((**iter).id()%1000==443||(**iter).id()%1000==553)&&
	 (**iter).children().size()==3) particles.push_back(*iter);
    }
  }
  // analyse them
  analyze(particles);
}

void OniumToOniumPiPiAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void OniumToOniumPiPiAnalysis::analyze(tPPtr part) {
  tPPtr out;
  Lorentz5Momentum q,pip;
  unsigned int npi0(0),npip(0),npim(0);
  for(unsigned int ix=0;ix<part->children().size();++ix) {
    if(part->children()[ix]->id()==ParticleID::piplus) {
      ++npip;
      pip=part->children()[ix]->momentum();
      q+=part->children()[ix]->momentum();
    }
    else if(part->children()[ix]->id()==ParticleID::piminus) {
      ++npim;
      q+=part->children()[ix]->momentum();
    }
    else if(part->children()[ix]->id()==ParticleID::pi0) {
      ++npi0;
      q+=part->children()[ix]->momentum();
      pip=part->children()[ix]->momentum();
    }
    else {
      out = part->children()[ix];
    }
  }
  // require pi+pi- or pi0pi0
  if(!(npi0==2||(npim==1&&npim==1))) return;
  // require incoming and outgoing onium
  if(part->id()%1000!=out->id()%1000) return;
  unsigned int ix=0; bool found(false);
  while(!found&&ix<_incoming.size()) {
    if(_incoming[ix]==part->id()&&_outgoing[ix]==out->id()) found=true;
    else ++ix;
  }
  if(!found) {
    Energy twompi=2.*getParticleData(ParticleID::piplus)->mass();
    Energy upp = part->nominalMass()+part->dataPtr()->widthUpCut()-
                 out ->nominalMass()+out->dataPtr()->widthLoCut();
    ix=_incoming.size();
    _incoming.push_back(part->id());
    _outgoing.push_back(out ->id());
    _mpipi.push_back(make_pair(new_ptr(Histogram(twompi/GeV,upp/GeV,200)),
			       new_ptr(Histogram(twompi/GeV,upp/GeV,200))));
    _hel  .push_back(make_pair(new_ptr(Histogram(-1.,1.,200)),
			       new_ptr(Histogram( 0.,1.,200))));
  }
  // calculate the mass
  q.rescaleMass();
  Energy mq=q.mass();
  // calculate the helicity angle
  Boost boost = -q.boostVector();
  Lorentz5Momentum qp=out->momentum();
  Lorentz5Momentum ppi=pip;
  qp.boost(boost);
  q.boost(boost);
  ppi.boost(boost);
  double cX=-ppi.vect().unit()*qp.vect().unit();
  if(npi0==2) {
    *_mpipi[ix].second+=mq/GeV;
    *_hel[ix].second+=abs(cX);
  }
  else {
    *_mpipi[ix].first+=mq/GeV;
    *_hel[ix].first+=cX;
  }
}

NoPIOClassDescription<OniumToOniumPiPiAnalysis> OniumToOniumPiPiAnalysis::initOniumToOniumPiPiAnalysis;
// Definition of the static class description member.

void OniumToOniumPiPiAnalysis::Init() {

  static ClassDocumentation<OniumToOniumPiPiAnalysis> documentation
    ("There is no documentation for the OniumToOniumPiPiAnalysis class");

}

void OniumToOniumPiPiAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    string title= 
      getParticleData(_incoming[ix])->PDGName() + " -> " +
      getParticleData(_outgoing[ix])->PDGName();
    string temp = "Mass distrubtion for pi+pi- in " + title;
    using namespace HistogramOptions;
    _mpipi[ix].first->topdrawOutput(output,Frame,"BLACK",temp);
    temp = "Mass distrubtion for pi0pi0 in " + title;
    _mpipi[ix].second->topdrawOutput(output,Frame,"BLACK",temp);
    temp = "Helicity angle for pi+pi- in " + title;
    _hel[ix].first->topdrawOutput(output,Frame,"BLACK",temp);
    temp = "Helicity angle for pi0pi0 in " + title;
    _hel[ix].second->topdrawOutput(output,Frame,"BLACK",temp);
  }
}
