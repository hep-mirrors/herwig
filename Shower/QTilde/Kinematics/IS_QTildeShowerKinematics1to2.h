// -*- C++ -*-
//
// IS_QTildeShowerKinematics1to2.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_IS_QTildeShowerKinematics1to2_H
#define HERWIG_IS_QTildeShowerKinematics1to2_H
//
// This is the declaration of the IS_QTildeShowerKinematics1to2 class.
//

#include "Herwig/Shower/QTilde/Kinematics/ShowerKinematics.h"
#include "Herwig/Shower/QTilde/SplittingFunctions/Sudakov1to2FormFactor.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  
 *  This (concrete) class provides the specific Intial State shower
 *  kinematics information.
 *
 *  @see ShowerKinematics
 *  @see FS_QTildeShowerKinematics1to2
 *  @see Decay_QTildeShowerKinematics1to2
 *  @see KinematicsReconstructor
 */
class IS_QTildeShowerKinematics1to2: public ShowerKinematics {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   *  Construct in terms of the basis states
   */
  IS_QTildeShowerKinematics1to2()= default;

  /**
   * The default constructor.
   */
  IS_QTildeShowerKinematics1to2(Energy scale, double z, double phi, Energy pt, tSudakovPtr sud);
  //@}

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
   * @param parent The branching particle
   * @param children The particles produced in the branching
   * @param partnerType The type of evolution partner
   */
  virtual void updateChildren( const tShowerParticlePtr parent, 
			       const ShowerParticleVector & children,
			       unsigned int pTscheme,
			       ShowerPartnerType partnerType) const;

  virtual void resetChildren( const tShowerParticlePtr parent, 
			      const ShowerParticleVector & children) const;

  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the 
   * KinematicsReconstructor.
   * @param parent The branching particle
   * @param children The particles produced in the branching
   * @param partnerType The type of evolution partner
   */
  virtual void updateParent( const tShowerParticlePtr parent, 
			     const ShowerParticleVector & children,
			     unsigned int pTscheme,
			     ShowerPartnerType partnerType) const;

  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the 
   * KinematicsReconstructor.
   */
  virtual void reconstructParent( const tShowerParticlePtr parent, 
				  const ParticleVector & children ) const;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   * @param theLast The particle.
   * @param px The \f$x\f$ component of the \f$p_T\f$.
   * @param py The \f$y\f$ component of the \f$p_T\f$.
   */
  virtual void updateLast(const tShowerParticlePtr theLast,
			  Energy px, Energy py) const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  IS_QTildeShowerKinematics1to2 & operator=(const IS_QTildeShowerKinematics1to2 &) = delete;

private:

  /**
   * \f$1\to2\f$ sudakov for extra features
   */
  tSudakov1to2Ptr sudakov1to2_;

};

}

#endif /* HERWIG_IS_QTildeShowerKinematics1to2_H */
