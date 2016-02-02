// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EvtGenDecayer class.
//

#include "EvtGenDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

bool EvtGenDecayer::accept(const DecayMode &) const {
  return true;
}

ParticleVector EvtGenDecayer::decay(const DecayMode & dm,
				    const Particle & parent) const {
  ParticleVector output;
  if(evtOpt_==0)
    output=evtgen_->decay(parent,false,dm);
  else if(evtOpt_==1)
    output=evtgen_->decay(parent, true,dm);
  else
    throw Exception() << "Unknown option in EvtGenDecayer::decay() " 
		      << Exception::runerror;
  if(check_) {
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

IBPtr EvtGenDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr EvtGenDecayer::fullclone() const {
  return new_ptr(*this);
}

void EvtGenDecayer::persistentOutput(PersistentOStream & os) const {
  os << evtgen_ << check_ << evtOpt_;
}

void EvtGenDecayer::persistentInput(PersistentIStream & is, int) {
  is >> evtgen_ >> check_ >> evtOpt_;  
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EvtGenDecayer,Decayer>
  describeThePEGEvtGenDecayer("Herwig::EvtGenDecayer", "HwEvtGenInterface.so");

void EvtGenDecayer::Init() {

  static ClassDocumentation<EvtGenDecayer> documentation
    ("The EvtGenDecayer class allows the EvtGen decay package to be used as"
     " a decayer inside Herwig");

  static Reference<EvtGenDecayer,EvtGenInterface> interfaceEvtGen
    ("EvtGen",
     "Pointer to the EvtGenInterface object which encapsulates the EvtGen decay package.",
     &EvtGenDecayer::evtgen_, false, false, true, false, false);

  static Switch<EvtGenDecayer,bool> interfaceCheck
    ("Check",
     "Perform some basic checks of the decay",
     &EvtGenDecayer::check_, false, false, false);
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

  static Switch<EvtGenDecayer,unsigned int> interfaceOption
    ("Option",
     "The way in which EvtGen is used.",
     &EvtGenDecayer::evtOpt_, 0, false, false);
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
