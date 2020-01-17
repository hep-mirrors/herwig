// -*- C++ -*-
//
// GeneralCurrentDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralCurrentDecayer class.
//

#include "GeneralCurrentDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void GeneralCurrentDecayer::persistentOutput(PersistentOStream & os) const {
  os << theVertex_ << inpart_ << outpart_ << currentOut_
     << current_ << ounit(maxmass_,GeV)
     << mode_ << wgtloc_ << wgtmax_ << weights_;
}

void GeneralCurrentDecayer::persistentInput(PersistentIStream & is, int) {
  is >> theVertex_ >> inpart_ >> outpart_ >> currentOut_
     >> current_ >> iunit(maxmass_,GeV)
     >> mode_ >> wgtloc_ >> wgtmax_ >> weights_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<GeneralCurrentDecayer,DecayIntegrator>
describeHerwigGeneralCurrentDecayer("Herwig::GeneralCurrentDecayer", "Herwig.so");

void GeneralCurrentDecayer::Init() {

  static ClassDocumentation<GeneralCurrentDecayer> documentation
    ("The GeneralCurrentDecayer class is designed to be the base class for all "
     "decays using the WeakDecayCurrents");
}

void GeneralCurrentDecayer::setDecayInfo(PDPtr in, PDPtr out, const vector<tPDPtr> & outCurrent,
					 VertexBasePtr vertex, WeakDecayCurrentPtr current,
					 Energy maxmass) {
  inpart_ = in;
  outpart_ = out;
  currentOut_ = outCurrent;
  theVertex_ = vertex;
  current_ = current;
  maxmass_ = maxmass;
}



int GeneralCurrentDecayer::modeNumber(bool & cc, tcPDPtr parent, 
				      const tPDVector & children) const {
  vector<long> id;
  id.push_back(parent->id());
  for(unsigned int ix=0;ix<children.size();++ix) id.push_back(children[ix]->id());
  return modeNumber(cc,id);
}

void GeneralCurrentDecayer::doinitrun() {
  current_->initrun();
  DecayIntegrator::doinitrun();
}

void GeneralCurrentDecayer::doinit() {
  DecayIntegrator::doinit();
  // make sure the current got initialised
  current_->init();
  // set up the phase-space channels
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr channel;
  tPDVector extpart;
  extpart.push_back(inpart_);
  extpart.push_back(outpart_);
  Energy mdiff(inpart_->mass()-outpart_->mass());
  vector<double> channelwgts;
  int iq(0),ia(0);
  bool done;
  vector<double>::iterator start,end;
  for(unsigned int ix=0;ix<current_->numberOfModes();++ix) {
    // get the external particles for this mode
    extpart.resize(2);
    tPDVector ptemp=current_->particles(inpart_->iCharge()-outpart_->iCharge(),ix,iq,ia);
    // check this is the right mode
    if(ptemp.size()!=currentOut_.size()) continue;
    vector<bool> matched(ptemp.size(),false);
    bool match = true;
    for(unsigned int iy=0;iy<currentOut_.size();++iy) {
      bool found = false;
      for(unsigned int iz=0;iz<ptemp.size();++iz) {
	if(!matched[iz]&&ptemp[iz]==currentOut_[iy]) {
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
    done=current_->createMode(inpart_->iCharge()-outpart_->iCharge(),
			      ix,mode,2,1,channel,mdiff);
    if(done) {
      // the maximum weight and the channel weights
      // the weights for the channel
      if(weights_.empty()) {
	weights_.resize(mode->numberChannels(),1./(mode->numberChannels()));
      }
      mode_ = ix;
      // special for the two body modes
      if(extpart.size()==3) {
	weights_.clear();
	mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
      }
      addMode(mode,wgtmax_,weights_);
    }
    break;
  }
}

int GeneralCurrentDecayer::modeNumber(bool & cc, vector<long> id) const {
  // incoming particle
  int idtemp[2];
  tPDPtr p0=getParticleData(id[0]);
  idtemp[0] = p0->CC() ? -id[0] : id[0];
  if(id    [0]==inpart_->id())      cc=false;
  else if(idtemp[0]==inpart_->id()) cc=true ;
  else return -1;
  tPDPtr p1 = outpart_;
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
  unsigned int icurr=current_->decayMode(idother);
  if(mode_==icurr) return  0;
  else             return -1;
}
