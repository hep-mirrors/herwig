// -*- C++ -*-
//
// Decay_QtildaShowerKinematics1to2.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Decay_QtildaShowerKinematics1to2_H
#define HERWIG_Decay_QtildaShowerKinematics1to2_H
//
// This is the declaration of the Decay_QtildaShowerKinematics1to2 class.
//

#include "QtildaShowerKinematics1to2.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This (concrete) class provides the specific decay shower
 *  kinematics information.
 *
 *  @see QtildaShowerKinematics1to2
 *  @see IS_QtildaShowerKinematics1to2
 *  @see FS_QtildaShowerKinematics1to2
 *  @see KinematicsReconstructor
 *
 */
class Decay_QtildaShowerKinematics1to2: public QtildaShowerKinematics1to2 {

public:

  /**
   *  The updateChildren, updateParent and updateLast
   *  members to update the values of the \f$\alpha\f$ and 
   *  \f$p_\perp\f$ variables during the shower evolution.
   */
  //@{

  /**
   * Along with the showering evolution --- going forward for
   * time-like (forward) evolution, and going backward for space-like
   * (backward) evolution --- the kinematical variables of the
   * branching products are calculated and updated from the knowledge
   * of the parent kinematics.  This method is used by the
   * ForwardShowerEvolver.  
   * ***ACHTUNG*** Might be extended to update colour connections as well.
   * @param theParent The branching particle
   * @param theChildren The particles produced in the branching
   */
  virtual void updateChildren( const tShowerParticlePtr theParent, 
			       const ShowerParticleVector & theChildren,
			       bool angularOrder ) const;

  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the 
   * KinematicsReconstructor.
   */
  virtual void reconstructParent( const tShowerParticlePtr theParent, 
				  const ParticleVector & theChildren ) const;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   * @param theLast The particle to update
   * @param iopt The option for the momentum reconstruction 
   * - 0 is in the rest frame of the pair of reference vectors
   * - 1 is in the rest frame of the p vector
   * @param mass The mass to be used, if less than zero on-shell
   */
  virtual void reconstructLast(const tShowerParticlePtr theLast,
			       unsigned int iopt, Energy mass=-1.*GeV) const;

  /**
   *  Perform any initial calculations needed after the branching has been selected
   * @param particle The branching particle
   * @param parent The bema particle for the jet if needed
   */
  virtual void initialize(ShowerParticle & particle,PPtr parent);
  //@}

private:
  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Decay_QtildaShowerKinematics1to2 & operator=(const Decay_QtildaShowerKinematics1to2 &);

};

}

#endif /* HERWIG_Decay_QtildaShowerKinematics1to2_H */
