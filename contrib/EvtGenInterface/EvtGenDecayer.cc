// -*- C++ -*-
//
// EvtGenDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EvtGenDecayer class.
//

#include "EvtGenDecayer.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/DecayMode.h"

namespace Herwig {
using namespace ThePEG;

bool EvtGenDecayer::accept(const DecayMode & dm) const {
  return true;
}

ParticleVector EvtGenDecayer::decay(const DecayMode & dm,
				    const Particle & parent) const {
  ParticleVector output;
  if(_evtopt==0)      output=_evtgen->decay(parent,false,dm);
  else if(_evtopt==1) output=_evtgen->decay(parent, true,dm);
  else  throw Exception() << "Unknown option in EvtGenDecayer::decay() " 
			  << Exception::runerror;
  if(_check) {
    Lorentz5Momentum ptotal=parent.momentum();
    int charge=parent.dataPtr()->iCharge();
    for(unsigned int ix=0;ix<output.size();++ix) {
      ptotal-=output[ix]->momentum();
      charge-=output[ix]->dataPtr()->iCharge();
      checkDecay(output[ix]);
    }
    if(abs(ptotal.x())>0.001*MeV||abs(ptotal.y())>0.001*MeV||
       abs(ptotal.z())>0.001*MeV||abs(ptotal.e())>0.001*MeV) {
      generator()->log() << "Decay of " << parent.PDGName() 
			 << " violates momentum conservation in"
			 << "EvtGenDecayer::decay\n";
    }
    if(charge!=0) {
      generator()->log() << "Decay of " << parent.PDGName() 
			 << " violates charge conservation in"
			 << "EvtGenDecayer::decay\n";
    }
  }
  return output;
}


void EvtGenDecayer::persistentOutput(PersistentOStream & os) const {
  os << _evtgen << _evtopt << _check;
}

void EvtGenDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _evtgen >> _evtopt >> _check;  
}

ClassDescription<EvtGenDecayer> EvtGenDecayer::initEvtGenDecayer;
// Definition of the static class description member.

void EvtGenDecayer::Init() {

  static ClassDocumentation<EvtGenDecayer> documentation
    ("The EvtGenDecayer class allows the EvtGen decay package to be used as"
     " a decayer inside Herwig++");

  static Reference<EvtGenDecayer,EvtGen> interfaceEvtGen
    ("EvtGen",
     "Pointer to the EvtGen object which encapsulates the EvtGen decay package.",
     &EvtGenDecayer::_evtgen, false, false, true, false, false);

  static Switch<EvtGenDecayer,unsigned int> interfaceOption
    ("Option",
     "The way in which EvtGen is used.",
     &EvtGenDecayer::_evtopt, 0, false, false);
  static SwitchOption interfaceOptionParent
    (interfaceOption,
     "Parent",
     "EvtGen decays the particle and returns the decay products to be decayed by"
     " Herwig++.",
     0);
  static SwitchOption interfaceOptionAll
    (interfaceOption,
     "All",
     "EvtGen decays the particle and all the unstable particles produced in the decay.",
     1);

  static Switch<EvtGenDecayer,bool> interfaceCheck
    ("Check",
     "Perform some basic checks of the decay",
     &EvtGenDecayer::_check, false, false, false);
  static SwitchOption interfaceCheckCheck
    (interfaceCheck,
     "Yes",
     "Perform the checks",
     true);
  static SwitchOption interfaceCheckNoCheck
    (interfaceCheck,
     "No",
     "Don't perform the checks",
     false);
}

void EvtGenDecayer::checkDecay(PPtr in) const {
  Lorentz5Momentum ptotal=in->momentum();
  int charge = in->dataPtr()->iCharge();
  if(in->children().empty()) return;
  for(unsigned int ix=0;ix<in->children().size();++ix) {
    checkDecay(in->children()[ix]);
    ptotal-=in->children()[ix]->momentum();
    charge-=in->children()[ix]->dataPtr()->iCharge();
  }
  if(abs(ptotal.x())>MeV||abs(ptotal.y())>MeV||
     abs(ptotal.z())>MeV||abs(ptotal.e())>MeV) {
    generator()->log() 
      << "Decay of " << in->PDGName() << " violates momentum conservation"
      << "in EvtGenDecayer::checkDecay";
  }
  if(charge!=0) generator()->log() << "Decay of " << in->PDGName() 
				   << " violates charge conservation in "
				   << "EvtGenDecayer::checkDecay\n";
}

}
