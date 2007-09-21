// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BranchingRatioAnalysis class.
//

#include "BranchingRatioAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

void BranchingRatioAnalysis::analyze(tEventPtr event, long , int , int ) {
  // find all vector particles 
  tPVector particles;
  for(unsigned int ix=0, nstep=event->primaryCollision()->steps().size();
      ix<nstep;++ix) {
    ThePEG::ParticleSet part=event->primaryCollision()->step(ix)->all();
    ThePEG::ParticleSet::iterator iter=part.begin();
    ThePEG::ParticleSet::iterator end=part.end();
    for( ;iter!=end;++iter) {
      if((**iter).id()==_particle->id()) particles.push_back(*iter);
    }
  }
  // analyse them
  analyze(particles);
}

void BranchingRatioAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void BranchingRatioAnalysis::analyze(tPPtr part) {
  tDMPtr mode=part->decayMode();
  if(!mode) return;
  if(_points.find(mode)==_points.end()) {
    _points[mode]=vector<double>(101,0.);
  }
  if(part->mass()>=part->dataPtr()->massMin()&&
     part->mass()<=part->dataPtr()->massMax()) {
    int ibin = 
      int((part->mass()              -part->dataPtr()->massMin())/
	  (part->dataPtr()->massMax()-part->dataPtr()->massMin())*100.);
    if(ibin>=0&&ibin<100) {
      _points[mode][ibin]+=1.;
      _total[ibin]+=1.;
      _points[mode][100]+=1.;
      _total[100]+=1.;
    }
  }
}

void BranchingRatioAnalysis::persistentOutput(PersistentOStream & os) const {
  os << _particle;
}

void BranchingRatioAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> _particle;
}

ClassDescription<BranchingRatioAnalysis> 
BranchingRatioAnalysis::initBranchingRatioAnalysis;
// Definition of the static class description member.

void BranchingRatioAnalysis::Init() {

  static ClassDocumentation<BranchingRatioAnalysis> documentation
    ("There is no documentation for the BranchingRatioAnalysis class");

  static Reference<BranchingRatioAnalysis,ParticleData> interfaceParticle
    ("Particle",
     "The particle for which the analysis is being perform",
     &BranchingRatioAnalysis::_particle, false, false, true, false, false);

}

void BranchingRatioAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + 
    string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  map<tDMPtr,vector<double> >::const_iterator it;
  Energy step = (_particle->massMax()-_particle->massMin())/100.;
  for(it=_points.begin();it!=_points.end();++it) {
    output << "NEW FRAME\n";
    output << "SET FONT DUPLEX\n";
    output << "TITLE TOP \"BR for " << it->first->tag() 
	   << "Total = " << it->second[100]/_total[100]  << "\"\n";
    output << "SET LIMITS X " 
	   << _particle->massMin()/GeV << " " 
	   << _particle->massMax()/GeV << "\n";
    Energy mass = _particle->massMin()-0.5*step;
    for(unsigned int ix=0;ix<it->second.size()-1;++ix) {
      mass+=step;
      if(_total[ix]>0.) {
	output << mass/GeV << "\t" << it->second[ix]/_total[ix] << "\n";
      }
      else {
	output << mass/GeV << "\t 0. \n";
      }
    }
    output << "HIST\n";
    output << "NEW FRAME\n";
    output << "SET FONT DUPLEX\n";
    output << "TITLE TOP \"Mass distribution for " 
	   << it->first->tag() << "\"\n";
    mass = _particle->massMin()-0.5*step;
    for(unsigned int ix=0;ix<it->second.size()-1;++ix) {
      mass+=step;
      if(_total[ix]>0.) {
	output << mass/GeV << "\t" 
	       << it->second[ix]/it->second[100]*GeV/step
	       << "\n";
      }
      else {
	output << mass/GeV << "\t 0. \n";
      }
    }
    output << "HIST\n";
  }
}
