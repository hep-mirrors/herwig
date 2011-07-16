// -*- C++ -*-
//
// GeneralCurrentDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
  os << _theVertex << _inpart << _outpart << _current << ounit(_maxmass,GeV)
     << _modemap << _modestart << _wgtloc 
     << _wgtmax << _weights;
}

void GeneralCurrentDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theVertex >> _inpart >> _outpart >> _current >> iunit(_maxmass,GeV)
     >> _modemap >> _modestart >> _wgtloc 
     >> _wgtmax >> _weights;
}

AbstractClassDescription<GeneralCurrentDecayer> GeneralCurrentDecayer::initGeneralCurrentDecayer;
// Definition of the static class description member.

void GeneralCurrentDecayer::Init() {

  static ClassDocumentation<GeneralCurrentDecayer> documentation
    ("The GeneralCurrentDecayer class is designed to be the base class for all "
     "decays using the WeakDecayCurrents");

  static Reference<GeneralCurrentDecayer,Helicity::VertexBase> interfaceDecayVertex
    ("DecayVertex",
     "Pointer to decayer vertex",
     &GeneralCurrentDecayer::_theVertex, false, false, true, false);

  static ParVector<GeneralCurrentDecayer,int> interfaceIncomingPart
    ("Incoming",
     "PDG Codes for incoming particles",
     &GeneralCurrentDecayer::_inpart, 0, -1, 0, 0,
     false, false, false);

  static ParVector<GeneralCurrentDecayer,int> interfaceOutgoingPart
    ("Outgoing",
     "PDG Codes for the outgoing particles",
     &GeneralCurrentDecayer::_outpart, 0, -1, 0, 0,
     false, false, false);

  static Reference<GeneralCurrentDecayer,WeakDecayCurrent> interfaceCurrent
    ("Current",
     "The weak current for the decay",
     &GeneralCurrentDecayer::_current, false, false, true, false, false);

  static Parameter<GeneralCurrentDecayer,Energy> interfaceMaximumMass
    ("MaximumMass",
     "The maximum mass difference for the decay",
     &GeneralCurrentDecayer::_maxmass, GeV, 5.0*GeV, 1.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
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
  _modemap.clear();
  _modestart.clear();
  // extract the possible particles for the modes
  vector<long> particles;
  for(unsigned int ix=0;ix<3;++ix) {
    vector<long>   temp      = _theVertex->search(ix,ParticleID::Wplus);
    particles.insert(particles.end(),temp.begin(),temp.end());
  }
  _inpart.clear();
  _outpart.clear();
  while(!particles.empty()) {
    vector<tPDPtr> part;
    for(unsigned int ix=0;ix<3;++ix) {
      if(abs(particles.back())!=ParticleID::Wplus) 
	part.push_back(getParticleData(particles.back()));
      particles.pop_back();
    }
    if(part[0]->mass()<part[1]->mass()) swap(part[0],part[1]);
    if(part[0]->CC()) part[0]=part[0]->CC();
    if(part[0]->mass()-part[1]->mass()>_maxmass) continue;
    // set up the phase-space channels
    DecayPhaseSpaceModePtr mode;
    DecayPhaseSpaceChannelPtr channel;
    tPDVector extpart,ptemp;
    extpart.push_back(part[0]);
    extpart.push_back(part[1]);
    Energy mdiff(part[0]->mass()-part[1]->mass());
    double maxweight;
    vector<double> channelwgts;
    int iq(0),ia(0);
    unsigned int ix,iy;
    bool done;
    vector<double>::iterator start,end;
    _modestart.push_back(_modemap.size());
    _inpart.push_back(part[0]->id());
    _outpart.push_back(part[1]->id());
    for(ix=0;ix<_current->numberOfModes();++ix) {
      // get the external particles for this mode
      extpart.resize(2);
      ptemp=_current->particles(part[0]->iCharge()-part[1]->iCharge(),ix,iq,ia);
      for(iy=0;iy<ptemp.size();++iy){extpart.push_back(ptemp[iy]);}
      // create the mode
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
      // create the first piece of the channel
      channel = new_ptr(DecayPhaseSpaceChannel(mode));
      channel->addIntermediate(extpart[0],0,0.0,-1,1);
      done=_current->createMode(part[0]->iCharge()-part[1]->iCharge(),
				ix,mode,2,1,channel,mdiff);
      if(done) {
	// the maximum weight and the channel weights
	// the maximum
	if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
	else                             maxweight=0.;
	// the weights for the channel
	if(_wgtloc.size()>numberModes()&&
	   _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
	  start=_weights.begin()+_wgtloc[numberModes()];
	  end  = start+mode->numberChannels();
	  channelwgts=vector<double>(start,end);
	}
	else {
	  channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
	}
	_modemap.push_back(ix);
	// special for the two body modes
	if(extpart.size()==3) {
	  channelwgts.clear();
	  mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
	}
	addMode(mode,maxweight,channelwgts);
      }
    }
  }
}

int GeneralCurrentDecayer::modeNumber(bool & cc, vector<long> id) const {
  // incoming particle
  int idtemp[2];
  tPDPtr p0=getParticleData(id[0]);
  idtemp[0] = p0->CC() ? -id[0] : id[0];
  // loop over possible particles
  for(unsigned int ipart=0;ipart<_inpart.size();++ipart) {
    if(id    [0]==_inpart[ipart]) {
      cc=false;
    }
    else if(idtemp[0]==_inpart[ipart]) {
      cc=true;
    }
    else
      continue;
    tPDPtr p1 = getParticleData(_outpart[ipart]);
    if(cc) p1 = p1->CC() ? p1->CC() : p1;
    // if this in the particles
    vector<long>::iterator iloc = std::find(++id.begin(), id.end(), p1->id());
    if(idtemp[0]==id[0]&&iloc==id.end()) {
      iloc = std::find(++id.begin(), id.end(), p1->CC()->id());
    }
    if(iloc==id.end()) continue;
    vector<int> idother;
    for(vector<long>::iterator it=++id.begin();it!=id.end();++it) {
      if(it!=iloc) idother.push_back(*it);
    }
    unsigned int icurr=_current->decayMode(idother);
    int imode,imax;
    if(ipart+1<_modestart.size()) imax=_modestart[ipart+1];
    else imax=_modemap.size();
    for(imode=_modestart[ipart];imode<imax;++imode) {
      if(_modemap[imode]==icurr) break;
    }
    if(imode>=imax) imode=-1;
    if(imode>=0) return imode;
  }
  return -1;
}
