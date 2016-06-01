// -*- C++ -*-
//
// RealEmissionProcess.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_RealEmissionProcess_H
#define HERWIG_RealEmissionProcess_H

#include "PerturbativeProcess.h"
#include "RealEmissionProcess.fh"

namespace Herwig {
using namespace ThePEG;

/**
 *  Simple struct for hard processes and decays
 */
class RealEmissionProcess : public PerturbativeProcess {

public:
  
  /**
   *  Constructor
   */
  RealEmissionProcess(PerturbativeProcessPtr born) {
    for(vector<pair<PPtr,tPerturbativeProcessPtr> >::iterator it=born->incoming().begin();
	it!=born->incoming().end();++it)
      bornIncoming().push_back(*it);
    for(vector<pair<PPtr, PerturbativeProcessPtr> >::iterator it=born->outgoing().begin();
	it!=born->outgoing().end();++it)
      bornOutgoing().push_back(*it);
  }
  
  /**
   *  The incoming particles
   */
  vector<pair<PPtr,tPerturbativeProcessPtr> > & bornIncoming() {
    return bornIncoming_;
  }
  
  /**
   *  The outgoing particles
   */
  vector<pair<PPtr, PerturbativeProcessPtr> > & bornOutgoing() {
    return bornOutgoing_;
  }
  
public:

  /**
   *  The emitter
   */
  unsigned int emitter  () const {return emitter_;}

  /**
   *  The spectator
   */
  unsigned int spectator() const {return spectator_;}

  /**
   *  The emitted
   */
  unsigned int emitted  () const {return emitted_;}

  /**
   *  The emitter
   */
  void emitter  (unsigned int in) {emitter_=in;}

  /**
   *  The spectator
   */
  void spectator(unsigned int in) {spectator_=in;}

  /**
   *  The emitted
   */
  void emitted  (unsigned int in) {emitted_=in;}

private:

  /**
   *  The emitter
   */
  unsigned int emitter_;

  /**
   *  The spectator
   */
  unsigned int spectator_;
  /**
   *  The emitter
   */
  unsigned int emitted_;

  /**
   *  Incoming for the Born process
   */
  vector<pair<PPtr,tPerturbativeProcessPtr> > bornIncoming_;

  /**
   *  Outgoing for the Born process
   */
  vector<pair<PPtr, PerturbativeProcessPtr> > bornOutgoing_;

};

}

#endif
