// -*- C++ -*-
#ifndef Herwig_FS_QTildeShowerKinematics1to1_H
#define Herwig_FS_QTildeShowerKinematics1to1_H
//
// This is the declaration of the FS_QTildeShowerKinematics1to1 class.
//

#include "ShowerKinematics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the FS_QTildeShowerKinematics1to1 class.
 */
class FS_QTildeShowerKinematics1to1: public ShowerKinematics {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FS_QTildeShowerKinematics1to1() = default;

  /**
   * The constructor.
   */
  FS_QTildeShowerKinematics1to1(Energy scale, tSudakovPtr sud) 
    : ShowerKinematics(scale,1.,0.,ZERO,sud) {}
  //@}

public:

  /**
   *
   */
  virtual void updateChildren( const tShowerParticlePtr parent, 
			       const ShowerParticleVector & children,
			       ShowerPartnerType partnerType) const;
  
  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   * @param parent   The parent
   * @param children The children
   * @param partnerType The type of evolution partner
   */
  virtual void updateParent(const tShowerParticlePtr parent,
			    const ShowerParticleVector & children,
			    unsigned int pTscheme,
			    ShowerPartnerType partnerType) const;

  virtual void resetChildren( const tShowerParticlePtr parent, 
			      const ShowerParticleVector & children) const;
public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FS_QTildeShowerKinematics1to1 & operator=(const FS_QTildeShowerKinematics1to1 &);

};

}

#endif /* Herwig_FS_QTildeShowerKinematics1to1_H */
