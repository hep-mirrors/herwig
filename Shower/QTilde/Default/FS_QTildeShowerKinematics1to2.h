// -*- C++ -*-
//
// FS_QTildeShowerKinematics1to2.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_FS_QTildeShowerKinematics1to2_H
#define HERWIG_FS_QTildeShowerKinematics1to2_H
//
// This is the declaration of the FS_QTildeShowerKinematics1to2 class.
//

#include "Herwig/Shower/Core/Base/ShowerKinematics.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This (concrete) class provides the specific Final State shower
 *  kinematics information.
 *
 *  @see ShowerKinematics
 *  @see IS_QTildeShowerKinematics1to2
 *  @see Decay_QTildeShowerKinematics1to2
 *  @see KinematicsReconstructor
 */
class FS_QTildeShowerKinematics1to2: public ShowerKinematics {

public:

  /**
   * Default constructor
   */
  inline FS_QTildeShowerKinematics1to2() {}

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
private:

  void updateParameters(tShowerParticlePtr theParent,
			tShowerParticlePtr theChild0,
			tShowerParticlePtr theChild1,
			bool setAlpha) const; 

public:
  virtual void updateChildren( const tShowerParticlePtr parent, 
			       const ShowerParticleVector & children,
			       ShowerPartnerType partnerType,
			       bool massVeto ) const;

  virtual void resetChildren( const tShowerParticlePtr parent, 
			      const ShowerParticleVector & children) const;


  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   * @param parent   The parent
   * @param children The children
   * @param partnerType The type of evolution partner
   */
  virtual void updateParent(const tShowerParticlePtr parent,
			    const ShowerParticleVector & children,
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
   * @param last The particle to update
   * @param mass The mass to be used, if less than zero on-shell
   */
  virtual void reconstructLast(const tShowerParticlePtr last, Energy mass=-1.*GeV) const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FS_QTildeShowerKinematics1to2 & operator=(const FS_QTildeShowerKinematics1to2 &);

};

}

#endif /* HERWIG_FS_QTildeShowerKinematics1to2_H */
