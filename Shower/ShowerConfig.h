// -*- C++ -*-
//
// ShowerConfig.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerConfig_H
#define HERWIG_ShowerConfig_H
//
// This is the declaration of the ShowerConfig class.

#include "Herwig++/Config/Herwig.h"
#include "Base/ShowerParticle.fh"
#include "Base/ShowerKinematics.fh" 

namespace Herwig { 
using namespace ThePEG;

/** \ingroup Shower
 *  
 *  Handy header file to be included in all Shower classes.
 *  It contains only some useful typedefs.
 */

  /**
   *  Forward declaration of the ColourLine of ThePEG
   */
  class ThePEG::ColourLine;
  
  /**
   * Pointer to a ColourLine
   */
  typedef Ptr<ThePEG::ColourLine>::pointer ColinePtr;
  
  /**
   * Transient Pointer to a ColourLine
   */
  typedef Ptr<ThePEG::ColourLine>::transient_pointer tColinePtr;
  
  /**
   * A pair of ColourLine pointers
   */
  typedef pair<ColinePtr,ColinePtr> ColinePair;

  /**
   * A pair of transient ColourLine pointers
   */
  typedef pair<tColinePtr,tColinePtr> tColinePair;

  /**
   *  A Vector of ShowerParticle pointers
   */
  typedef vector<ShowerParticlePtr> ShowerParticleVector;

  /**
   *  A Vector of transient ShowerParticle pointers
   */
  typedef vector<tShowerParticlePtr> tShowerParticleVector;

  /**
   *  Definition of the IdList for branchings
   */
  typedef vector<long> IdList;

}

#endif // HERWIG_ShowerConfig_H 

