// -*- C++ -*-
//
// DipolePartonSplitter.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipolePartonSplitter_H
#define HERWIG_DipolePartonSplitter_H
//
// This is the declaration of the DipolePartonSplitter class.
//

#include "ThePEG/EventRecord/Particle.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 * 
 * \brief The DipolePartonSplitter is a helper class to fix up
 * colour and mother-child relations in typical shower
 * splittings.
 *
 */
struct DipolePartonSplitter {

  /**
   * Fix up mother child relations for splitting the first
   * to the second and third parton; use conventions
   * of a backward splitting, if initialState is true.
   * In this case, the first child is assumed to be the
   * new incoming parton.
   */
  static void split(tPPtr parent, tPPtr firstChild, tPPtr secondChild, 
		    bool initialState, bool decayedEmitter=false);

  /**
   * Fix up relations for splitting the first
   * to the second and third parton; use conventions
   * of a backward splitting, if initialState is true.
   * In this case, the first child is assumed to be the
   * new incoming parton. A reference is given to
   * determine which colour line did actually radiate
   * in case of ambiguities, in the sense that the
   * colour connected parent-ref pair raidated as a dipole.
   */
  static void split(tPPtr parent, tPPtr firstChild, tPPtr secondChild, 
		    tPPtr ref, bool initialState, bool decayedEmitter=false);

  /**
   * Fix up relations for the case that
   * the new parton instance exists only
   * due to changes in e.g. kinematics.
   */
  static void change(tPPtr parent, tPPtr child, bool initialState, bool decayedSpec=false);

  /**
   * Return true, if the given partons are colour connected.
   */
  static bool colourConnected(tcPPtr first, tcPPtr second);

};

}

#endif /* HERWIG_DipolePartonSplitter_H */
