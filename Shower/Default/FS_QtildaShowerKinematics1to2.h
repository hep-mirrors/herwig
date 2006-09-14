// -*- C++ -*-
#ifndef HERWIG_FS_QtildaShowerKinematics1to2_H
#define HERWIG_FS_QtildaShowerKinematics1to2_H
//
// This is the declaration of the FS_QtildaShowerKinematics1to2 class.
//

#include "QtildaShowerKinematics1to2.h"
#include "FS_QtildaShowerKinematics1to2.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This (concrete) class provides the specific Final State shower
 *  kinematics information.
 *
 *  @see QtildaShowerKinematics1to2
 *  @see IS_QtildaShowerKinematics1to2
 *  @see KinematicsReconstructor
 */
class FS_QtildaShowerKinematics1to2: public QtildaShowerKinematics1to2 {

public:

  /**
   * Default constructor
   */
  inline FS_QtildaShowerKinematics1to2();

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
			       const ShowerParticleVector theChildren ) const;

  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the 
   * KinematicsReconstructor.
   */
  virtual void updateParent( const tShowerParticlePtr theParent, 
			     const ParticleVector theChildren ) const;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   * @param theLast The particle to update
   * @param iopt The option for the momentum reconstruction 
   * - 0 is in the rest frame of the pair of reference vectors
   * - 1 is in the rest frame of the p vector
   */
  virtual void updateLast( const tShowerParticlePtr theLast,unsigned int iopt ) const;

  /**
   *  Perform any initial calculations needed after the branching has been selected
   */
  virtual void initialize(ShowerParticle & particle);
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FS_QtildaShowerKinematics1to2 & operator=(const FS_QtildaShowerKinematics1to2 &);

};

}

#include "FS_QtildaShowerKinematics1to2.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FS_QtildaShowerKinematics1to2.tcc"
#endif

#endif /* HERWIG_FS_QtildaShowerKinematics1to2_H */
