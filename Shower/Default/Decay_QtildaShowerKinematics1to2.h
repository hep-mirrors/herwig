// -*- C++ -*-
#ifndef HERWIG_Decay_QtildaShowerKinematics1to2_H
#define HERWIG_Decay_QtildaShowerKinematics1to2_H
//
// This is the declaration of the Decay_QtildaShowerKinematics1to2 class.
//

#include "QtildaShowerKinematics1to2.h"
#include "Decay_QtildaShowerKinematics1to2.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the Decay_QtildaShowerKinematics1to2 class.
 */
class Decay_QtildaShowerKinematics1to2: public QtildaShowerKinematics1to2 {

public:

  /**
   * Default constructor 
   */
  inline Decay_QtildaShowerKinematics1to2();

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
  virtual void reconstructParent( const tShowerParticlePtr theParent, 
				  const ParticleVector theChildren ) const;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   */
  virtual void reconstructLast( const tShowerParticlePtr theLast,unsigned int iopt) const;

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
  Decay_QtildaShowerKinematics1to2 & operator=(const Decay_QtildaShowerKinematics1to2 &);

};

}

#include "Decay_QtildaShowerKinematics1to2.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Decay_QtildaShowerKinematics1to2.tcc"
#endif

#endif /* HERWIG_Decay_QtildaShowerKinematics1to2_H */
