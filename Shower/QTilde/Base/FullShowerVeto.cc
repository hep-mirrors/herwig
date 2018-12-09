// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FullShowerVeto class.
//

#include "FullShowerVeto.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void FullShowerVeto::persistentOutput(PersistentOStream & os) const {
  os << type_ << behaviour_;
}

void FullShowerVeto::persistentInput(PersistentIStream & is, int) {
  is >> type_ >> behaviour_;
}

void FullShowerVeto::doinit() {
  Interfaced::doinit();
  // check reweighting
  if(behaviour_==1&&type_!=2)
    throw Exception() << "Reweighting in the FullShowerVeto is only"
		      << " supported for the primary hard process\n"
		      << Exception::runerror;
}

DescribeAbstractClass<FullShowerVeto,Interfaced>
describeHerwigFullShowerVeto("Herwig::FullShowerVeto", "HwShower.so");

void FullShowerVeto::Init() {

  static ClassDocumentation<FullShowerVeto> documentation
    ("The FullShowerVeto class allows the parton shower generated from a configuration to be vetoed.");

  static Switch<FullShowerVeto,unsigned int> interfaceType
    ("Type",
     "Which type of processes to consider",
     &FullShowerVeto::type_, 1, false, false);
  static SwitchOption interfaceTypeAll
    (interfaceType,
     "All",
     "All Processes",
     0);
  static SwitchOption interfaceTypeScattering
    (interfaceType,
     "Scattering",
     "Only apply to scattering processes and not decays",
     1);
  static SwitchOption interfaceTypePrimary
    (interfaceType,
     "Primary",
     "Only apply to the primary scattering process",
     2);
  static SwitchOption interfaceTypeDecay
    (interfaceType,
     "Decay",
     "Only apply to decays",
     3);

  static Switch<FullShowerVeto,unsigned int> interfaceBehaviour
    ("Behaviour",
     "What to do if the shower if vetoed",
     &FullShowerVeto::behaviour_, 0, false, false);
  static SwitchOption interfaceBehaviourShower
    (interfaceBehaviour,
     "Shower",
     "Veto the shower and try showering the process again",
     0);
  static SwitchOption interfaceBehaviourShowerReweight
    (interfaceBehaviour,
     "ShowerReweight",
     "Veto the shower and reweight the event to take this into account, only supported for the primary process",
     1);
  static SwitchOption interfaceBehaviourEvent
    (interfaceBehaviour,
     "Event",
     "Veto the event, cross section automatically reweigted",
     2);

}

int FullShowerVeto::applyVeto(ShowerTreePtr tree) {
  // return if veto should not be calculated
  // decay process and only doing hard processes
  // or vice versa
  if(((type_ == 1 || type_ ==2 ) && tree->isDecay() ) ||
     ( type_ == 3                && tree->isHard()))
    return -1;
  // not primary process and only doing those
  if( type_ == 2 && !ShowerHandler::currentHandler()->firstInteraction() )
     return -1;
  // extract the incoming and outgoing particles from the ShowerTree
  finalState_.clear();
  incoming_.clear();
  outgoing_.clear();
  // incoming
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator it=tree->incomingLines().begin();
      it!=tree->incomingLines().end();++it) {
    incoming_.push_back(it->first->progenitor());
  }
  // outgoing
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator it=tree->outgoingLines().begin();
      it!=tree->outgoingLines().end();++it) {
    outgoing_.push_back(it->first->progenitor());
  }
  // call function in inheriting class to decide what to do
  bool vetoed = vetoShower();
  // clear storage
  finalState_.clear();
  incoming_.clear();
  outgoing_.clear();
  // return the answer
  return vetoed ? int(behaviour_) : -1;
}

namespace {

void addFinal(vector<tPPtr> & finalState, tPPtr particle) {
  if(particle->children().empty()) {
    finalState.push_back(particle);
    return;
  }
  for(unsigned int ix=0;ix<particle->children().size();++ix) {
    addFinal(finalState,particle->children()[ix]);
  }
}
}


// extract the incoming and outgoing particles
const vector<tPPtr> & FullShowerVeto::finalState() {
  if(!finalState_.empty()) return finalState_;
  // incoming
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    if(incoming_[ix]->parents().empty()) continue;
    tPPtr parent = incoming_[ix]->parents()[0];
    while (parent) {
      addFinal(finalState_,parent->children()[1]);
      if(parent->parents().empty()) break;
      parent = parent->parents()[0];
    }
  }
  // outgoing
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    addFinal(finalState_,outgoing_[ix]);
  } 
  return finalState_;
}
