// -*- C++ -*-
//
// GeneralCurrentDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralCurrentDecayer class.
//

#include "GeneralCurrentDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void GeneralCurrentDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theVertex << _inpart << _outpart << _currentOut
     << _current << ounit(_maxmass,GeV)
     << _mode << _wgtloc << _wgtmax << _weights;
}

void GeneralCurrentDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theVertex >> _inpart >> _outpart >> _currentOut
     >> _current >> iunit(_maxmass,GeV)
     >> _mode >> _wgtloc >> _wgtmax >> _weights;
}

AbstractClassDescription<GeneralCurrentDecayer> GeneralCurrentDecayer::initGeneralCurrentDecayer;
// Definition of the static class description member.

void GeneralCurrentDecayer::Init() {

  static ClassDocumentation<GeneralCurrentDecayer> documentation
    ("The GeneralCurrentDecayer class is designed to be the base class for all "
     "decays using the WeakDecayCurrents");
}

void GeneralCurrentDecayer::setDecayInfo(PDPtr in, PDPtr out, const vector<tPDPtr> & outCurrent,
					 VertexBasePtr vertex, WeakDecayCurrentPtr current,
					 Energy maxmass) {
  _inpart = in;
  _outpart = out;
  _currentOut = outCurrent;
  _theVertex = vertex;
  _current = current;
  _maxmass = maxmass;
}



int GeneralCurrentDecayer::modeNumber(bool & cc, tcPDPtr parent, 
				      const tPDVector & children) const {
  vector<long> id;
  id.push_back(parent->id());
  for(unsigned int ix=0;ix<children.size();++ix) id.push_back(children[ix]->id());
  return modeNumber(cc,id);
}

void GeneralCurrentDecayer::doinitrun() {
  _current->initrun();
  DecayIntegrator::doinitrun();
}

void GeneralCurrentDecayer::doinit() {
  DecayIntegrator::doinit();
  // make sure the current got initialised
  _current->init();
  // set up the phase-space channels
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr channel;
  tPDVector extpart;
  extpart.push_back(_inpart);
  extpart.push_back(_outpart);
  Energy mdiff(_inpart->mass()-_outpart->mass());
  vector<double> channelwgts;
  int iq(0),ia(0);
  bool done;
  vector<double>::iterator start,end;
  for(unsigned int ix=0;ix<_current->numberOfModes();++ix) {
    // get the external particles for this mode
    extpart.resize(2);
    tPDVector ptemp=_current->particles(_inpart->iCharge()-_outpart->iCharge(),ix,iq,ia);
    // check this is the right mode
    if(ptemp.size()!=_currentOut.size()) continue;
    vector<bool> matched(ptemp.size(),false);
    bool match = true;
    for(unsigned int iy=0;iy<_currentOut.size();++iy) {
      bool found = false;
      for(unsigned int iz=0;iz<ptemp.size();++iz) {
	if(!matched[iz]&&ptemp[iz]==_currentOut[iy]) {
	  found = true;
	  matched[iz] = true;
	  break;
	}
      }
      if(!found) {
	match = false;
	break;
      }
    }
    if(!match) continue;
    for(unsigned int iy=0;iy<ptemp.size();++iy)
      extpart.push_back(ptemp[iy]);
    // create the mode
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    // create the first piece of the channel
    channel = new_ptr(DecayPhaseSpaceChannel(mode));
    channel->addIntermediate(extpart[0],0,0.0,-1,1);
    done=_current->createMode(_inpart->iCharge()-_outpart->iCharge(),
			      ix,mode,2,1,channel,mdiff);
    if(done) {
      // the maximum weight and the channel weights
      // the weights for the channel
      if(_weights.empty()) {
	_weights.resize(mode->numberChannels(),1./(mode->numberChannels()));
      }
      _mode = ix;
      // special for the two body modes
      if(extpart.size()==3) {
	_weights.clear();
	mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
      }
      addMode(mode,_wgtmax,_weights);
    }
    break;
  }
}

int GeneralCurrentDecayer::modeNumber(bool & cc, vector<long> id) const {
  // incoming particle
  int idtemp[2];
  tPDPtr p0=getParticleData(id[0]);
  idtemp[0] = p0->CC() ? -id[0] : id[0];
  if(id    [0]==_inpart->id())      cc=false;
  else if(idtemp[0]==_inpart->id()) cc=true ;
  else return -1;
  tPDPtr p1 = _outpart;
  if(cc&&p1->CC()) p1 = p1->CC();
  // if this in the particles
  vector<long>::iterator iloc = std::find(++id.begin(), id.end(), p1->id());
  if(idtemp[0]==id[0]&&iloc==id.end()) {
    iloc = std::find(++id.begin(), id.end(), p1->CC()->id());
  }
  if(iloc==id.end()) return -1;
  vector<int> idother;
  for(vector<long>::iterator it=++id.begin();it!=id.end();++it) {
    if(it!=iloc) idother.push_back(*it);
  }
  unsigned int icurr=_current->decayMode(idother);
  if(_mode==icurr) return  0;
  else             return -1;
}
