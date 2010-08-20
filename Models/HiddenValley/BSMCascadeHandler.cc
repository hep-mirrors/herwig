// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BSMCascadeHandler class.
//

#include "BSMCascadeHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/PDF/HwRemDecayer.h"
#include "AdditionalGaugeParticleData.h"

using namespace Herwig;

BSMCascadeHandler::BSMCascadeHandler() {}

BSMCascadeHandler::~BSMCascadeHandler() {}

IBPtr BSMCascadeHandler::clone() const {
  return new_ptr(*this);
}

IBPtr BSMCascadeHandler::fullclone() const {
  return new_ptr(*this);
}

void BSMCascadeHandler::persistentOutput(PersistentOStream & ) const {
}

void BSMCascadeHandler::persistentInput(PersistentIStream & , int) {
}

ClassDescription<BSMCascadeHandler> BSMCascadeHandler::initBSMCascadeHandler;
// Definition of the static class description member.

void BSMCascadeHandler::Init() {

  static ClassDocumentation<BSMCascadeHandler> documentation
    ("There is no documentation for the BSMCascadeHandler class");

}

namespace {

  void findChildren(PPtr parent, set<PPtr> & inter, set<PPtr> & final) {
    for(unsigned int ix=0;ix<parent->children().size();++ix) {
      if(parent->children()[ix]->children().empty())
	final.insert(parent->children()[ix]);
      else {
	inter.insert(parent->children()[ix]);
	findChildren(parent->children()[ix],inter,final);
      }
    }
  }
}


void BSMCascadeHandler::cascade() {
  setCurrentHandler();
  map<PPtr,ParticleVector> decaymap;
  for(unsigned int ix=0;ix<tagged().size();++ix) {
    tcAdditionalGaugeParticleDataPtr ptemp 
      = dynamic_ptr_cast<tcAdditionalGaugeParticleDataPtr>(tagged()[ix]->dataPtr());
    if(!ptemp) continue;
    PPtr parent = tagged()[ix];
    do {
      parent = parent->parents()[0];
      ptemp = dynamic_ptr_cast<tcAdditionalGaugeParticleDataPtr>(parent->dataPtr());
    }
    while(ptemp);
    map<PPtr,ParticleVector>::iterator it = decaymap.find(parent);
    if(it!=decaymap.end()) it->second.push_back(tagged()[ix]);
    else decaymap.insert(make_pair(parent,ParticleVector(1,tagged()[ix])));
  }
  for(map<PPtr,ParticleVector>::iterator it=decaymap.begin();
      it!=decaymap.end();++it) {
    vector<pair<tPPtr,tcAdditionalGaugeParticleDataPtr> > darkParticles;
    for(unsigned int ix=0;ix<it->second.size();++ix) {
      tcAdditionalGaugeParticleDataPtr ptemp 
	= dynamic_ptr_cast<tcAdditionalGaugeParticleDataPtr>(it->second[ix]->dataPtr());
      if(ptemp) darkParticles.push_back(make_pair(it->second[ix],ptemp));
    }
    if(darkParticles.size()==2) {
      ColinePtr newline(new_ptr(ColourLine()));
      for(unsigned int ix=0;ix<darkParticles.size();++ix) {
	if(darkParticles[ix].second->hiddenColour()==
	   HiddenPDT::HiddenColourFundamental)
	  newline->addColoured(darkParticles[ix].first);
	else if(darkParticles[ix].second->hiddenColour()==
		HiddenPDT::HiddenColourAntiFundamental)
	  newline->addAntiColoured(darkParticles[ix].first);
      }
    }
    else {
      cerr << "can't sort out the dark colour\n";
      exit(0);
    }
    ParticleVector children=it->first->children();
    for(unsigned int ix=0;ix<children.size();++ix) {
      it->first->abandonChild(children[ix]);
    }
    for(unsigned int ix=0;ix<it->second.size();++ix) {
      tParticleVector parents = it->second[ix]->parents();
      for(unsigned int iy=0;iy<parents.size();++iy)
	parents[iy]->abandonChild(it->second[ix]);
      it->first->addChild(it->second[ix]);
    }
    ShowerDecayMap decays;
    ShowerTreePtr newtree(new_ptr(ShowerTree(it->first,decays)));
    Energy width=it->first->dataPtr()->generateWidth(it->first->mass());
    decays.insert(make_pair(width,newtree));
    vector<ShowerTreePtr> done;
    while(!decays.empty()) {
      // find particle whose production process has been showered
      ShowerDecayMap::iterator dit = decays.begin();
      //while(!dit->second->parent()->hasShowered() && dit!=decays.end()) ++dit;
      assert(dit!=decays.end());
      // get the particle
      ShowerTreePtr decayingTree = dit->second;
      // remove it from the multimap
      decays.erase(dit);
      // make sure the particle has been decayed
      //decayingTree->decay(decays);
      // now shower the decay
      evolver()->showerDecay(decayingTree);
      done.push_back(decayingTree);
      decayingTree->updateAfterShower(decays);
    }
    for(unsigned int ix=0;ix<done.size();++ix) {
      done[ix]->fillEventRecord(currentStep(),false,true);
    }
    PPtr parent=it->first;
    while(parent->children()[0]->id()==parent->id()) {
      parent = parent->children()[0];
    }
    for(unsigned int ix=0;ix<children.size();++ix) {
      ParticleVector temp = children[ix]->children();
      for(unsigned int iy=0;iy<temp.size();++iy)
	currentStep()->removeParticle(temp[iy]);
    }
    ParticleVector newchildren=parent->children();
    for(unsigned int ix=0;ix<children.size();++ix)
      currentStep()->insertIntermediate(children[ix],parent,newchildren[ix]);
  }
}
