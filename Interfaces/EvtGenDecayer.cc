// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EvtGenDecayer class.
//

#include "EvtGenDecayer.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EvtGenDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"

namespace Herwig {
using namespace ThePEG;

EvtGenDecayer::~EvtGenDecayer() {}

bool EvtGenDecayer::accept(const DecayMode & dm) const {
  // this decayer should only be used if letting EvtGen pick the mode
  if((_evtopt==0||_evtopt==2)&&dm.wildProductMatcher()&&
     dm.wildProductMatcher()->name()=="MatchAny")
    {return true;}
  else if(_evtopt==1){return true;}
  return false;
}

ParticleVector EvtGenDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector output;
  if(_evtopt==0){output=_evtgen->randomDecayAll(parent);}
  else if(_evtopt==1){output=_evtgen->decayAll(dm,parent);}
  else if(_evtopt==2){output=_evtgen->randomDecay(parent);}
  else {throw Exception() << "Unknown option in EvtGenDecayer::decay() " 
			  << Exception::runerror;}
  return output;
}


void EvtGenDecayer::persistentOutput(PersistentOStream & os) const {
  os << _evtgen << _evtopt;
}

void EvtGenDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _evtgen >> _evtopt;  
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
  static SwitchOption interfaceOptionEvtGenAll
    (interfaceOption,
     "EvtGenAll",
     "vtGen selects the decay mode of the particle and decays all"
     " the unstable particles produced in the decay.",
     0);
  static SwitchOption interfaceOptionHerwigModeEvtGenAll
    (interfaceOption,
     "HerwigModeEvtGenAll",
     "Herwig++ selects the decay mode and then EvtGen decays the particle"
     " and all the unstable particles produced in the decay.",
     1);
  static SwitchOption interfaceOptionEvtGenMode
    (interfaceOption,
     "EvtGenMode",
     "EvtGen selects the decay mode and then EvtGen decays the particle"
     "but not the unstable particles produced in the decay.",
     2);
}

}
