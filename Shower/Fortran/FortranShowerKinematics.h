// -*- C++ -*-
#ifndef HERWIG_FortranShowerKinematics_H
#define HERWIG_FortranShowerKinematics_H
//
// This is the declaration of the FortranShowerKinematics class.
//

#include "Herwig++/Shower/Base/ShowerKinematics.h"
#include "FortranShowerKinematics.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the FortranShowerKinematics class.
 */
class FortranShowerKinematics: public ShowerKinematics {

public:

  /**
   * The default constructor.
   */
  inline FortranShowerKinematics();

  /**
   *  The updateChildren and updateParent
   *  members to update the values of the \f$\alpha\f$ and 
   *  \f$p_\perp\f$ variables during the shower evolution.
   */
  //@{
  /**
   * Along with the showering evolution --- going forward for
   * time-like (forward) evolution, and going backward for space-like
   * (backward) evolution --- the kinematical variables of the
   * branching products are calculated and updated from the knowledge
   * of the parent kinematics. 
   * @param theParent   The parent
   * @param theChildren The children
   */
  virtual void updateChildren(const tShowerParticlePtr theParent, 
			      const ShowerParticleVector theChildren,
			      bool angularOrder) const;

  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   * @param theParent   The parent
   * @param theChildren The children
   */
  virtual void updateParent(const tShowerParticlePtr theParent, 
			    const ShowerParticleVector theChildren) const;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   * @param theLast The particle.
   */
  virtual void updateLast(const tShowerParticlePtr theLast) const;
  //@}

  /**
   *  The reconstructLast, reconstructChildren and reconstructParent members
   *  are used during the reconstruction 
   */
  //@{
  /**
   * Along with the showering evolution --- going forward for
   * time-like (forward) evolution, and going backward for space-like
   * (backward) evolution --- the kinematical variables of the
   * branching products are calculated and updated from the knowledge
   * of the parent kinematics. 
   * @param theParent   The parent
   * @param theChildren The children
   */
  virtual void reconstructChildren(const tShowerParticlePtr theParent, 
			      const ShowerParticleVector theChildren) const;

  /**
   * Reconstruct the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   * @param theParent   The parent
   * @param theChildren The children
   */
  virtual void reconstructParent(const tShowerParticlePtr theParent, 
				 const ParticleVector theChildren) const;

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   * @param theLast The particle.
   * @param iopt The option for the momentum reconstruction 
   * - 0 is in the rest frame of the pair of reference vectors
   * - 1 is in the rest frame of the p vector
   */
  virtual void reconstructLast(const tShowerParticlePtr theLast,
			       unsigned int iopt) const;

  /**
   *  Perform any initial calculations needed after the branching has been selected
   */
  virtual void initialize(ShowerParticle & particle);
  //@}

  /**
   *  Don't need the basis vectors, throws exception
   */
  virtual vector<Lorentz5Momentum> getBasis() const; 

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FortranShowerKinematics & operator=(const FortranShowerKinematics &);

};

}

#include "FortranShowerKinematics.icc"

#endif /* HERWIG_FortranShowerKinematics_H */
