// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhaseSpaceChannel class.
//

#include "PhaseSpaceChannel.h"
#include "PhaseSpaceMode.h"

void PhaseSpaceChannel::init(tPhaseSpaceModePtr mode) {
  mode_=mode;
  // find the descentents
    for(PhaseSpaceResonance & res : intermediates_)
      findChildren(res,res.descendents);
    // ensure intermediates either have the width set, or
    // can't possibly be on-shell
    // first the maximum energy release
    Energy massmax = mode->eMax_;
    for(tcPDPtr part : mode->outgoing_)
      massmax -= mode->testOnShell_ ? part->mass() : part->massMin();
    for(PhaseSpaceResonance & res : intermediates_) {
      if(!res.particle || res.mWidth!=ZERO ||
	 res.jacobian != PhaseSpaceResonance::BreitWigner) continue;
      Energy massmin(ZERO);
      for(const int & ext : res.descendents) 
 	massmin += mode->testOnShell_ ? 
 	  mode->outgoing_[ext]->mass() : mode->outgoing_[ext]->massMin();
      // check if can be on-shell
      Energy mass = sqrt(res.mass2);
      if(mass>=massmin&&mass<=massmax+massmin) {
	string modeout = mode->incoming_.first->PDGName() + " ";
	if(mode->incoming_.second)
	  modeout += mode->incoming_.second->PDGName() + " ";
	for( tcPDPtr out : mode->outgoing_)
 	  modeout += out->PDGName() + " ";
 	throw InitException() << "Width zero for " << res.particle->PDGName()
 			      << " in PhaseSpaceChannel::init() "
 			      << modeout << Exception::runerror;
      }
    }
  }
