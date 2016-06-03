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

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Config/Pointers.h"
#include "RealEmissionProcess.fh"
#include "Herwig/Shower/QTilde/ShowerInteraction.h"

namespace Herwig {
using namespace ThePEG;

/**
 *  Simple struct for hard processes and decays
 */
class RealEmissionProcess : public Base {

public:

  /**
   *  The incoming particles
   */
  ParticleVector & incoming() {
    return incoming_;
  }

  /**
   *  The outgoing particles
   */
  ParticleVector & outgoing() {
    return outgoing_;
  }

public:
  
  /**
   *  The incoming particles
   */
  ParticleVector & bornIncoming() {
    return bornIncoming_;
  }
  
  /**
   *  The outgoing particles
   */
  ParticleVector & bornOutgoing() {
    return bornOutgoing_;
  }
  
  /**
   *  The hadrons
   */
  ParticleVector & hadrons() {
    return hadrons_;
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

public:

  /**
   *  Lorentz Rotation to final-state for II dipoles
   */ 
  LorentzRotation transformation() const { return trans_;}

  /**
   *  Lorentz Rotation to final-state for II dipoles
   */ 
  void transformation(LorentzRotation in) {trans_=in;}

public:

  /**
   *  Get the x values
   */
  pair<double,double> x() const {return x_;}
    
  /**
   *  Set the x values
   */
  void x(pair<double,double> in) {x_=in;}

public:

  /**
   *  Type of interaction
   */
  ShowerInteraction::Type interaction() {return interaction_;}

  /**
   *  Type of interaction
   */
  void interaction(ShowerInteraction::Type in) {interaction_ = in;}

  /**
   *  Emission scales
   */
  map<ShowerInteraction::Type,Energy> & pT() {return pT_;}

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

private:

  /**
   *  The incoming particles
   */
  ParticleVector incoming_;

  /*
   *   The outgoing particles
   */
  ParticleVector outgoing_;

private:

  /**
   *  The hadrons
   */
  ParticleVector hadrons_;

private:

  /**
   *  Incoming for the Born process
   */
  ParticleVector bornIncoming_;

  /**
   *  Outgoing for the Born process
   */
  ParticleVector bornOutgoing_;

private:

  /**
   *  Lorentz transformation for spectators in II
   */
  LorentzRotation trans_;

  /**
   *  x values
   */
  pair<double,double> x_;

private:

  /**
   *  Type of interaction
   */
  ShowerInteraction::Type interaction_;

  /**
   *  Emission scales
   */
  map<ShowerInteraction::Type,Energy> pT_;
};

}

#endif
