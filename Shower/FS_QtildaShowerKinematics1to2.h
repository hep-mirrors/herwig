// -*- C++ -*-
#ifndef HERWIG_FS_QtildaShowerKinematics1to2_H
#define HERWIG_FS_QtildaShowerKinematics1to2_H
//
// This is the declaration of the FS_QtildaShowerKinematics1to2 class.

#include "QtildaShowerKinematics1to2.h"


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
class FS_QtildaShowerKinematics1to2 : public QtildaShowerKinematics1to2 {

public:

  /**
   * Standard ctors and dtor.
   */
  inline FS_QtildaShowerKinematics1to2();
  inline FS_QtildaShowerKinematics1to2(const FS_QtildaShowerKinematics1to2 &);
  virtual ~FS_QtildaShowerKinematics1to2();

  /**
   * Creator with the two defining vectors  p  and  n . 
   */
  inline FS_QtildaShowerKinematics1to2( const Lorentz5Momentum & p, 
					const Lorentz5Momentum & n );

  /**
   * Along with the showering evolution --- going forward for
   * time-like (forward) evolution, and going backward for space-like
   * (backward) evolution --- the Sudakov variables associated to the
   * branching products are calcalted and returned, from the knowledge
   * of the parent Sudakov variables.   
   * Note that only alpha and p_perp can be reconstructed 
   * at this moment and we will obtain instead beta only later, 
   * using updateParent().
   */
  virtual void updateChildren( const double parentSudAlpha, 
			       const Energy parentSudPx, const Energy parentSudPy, 
                               vector<double> & sudAlphaVect, 
			       vector<Energy> & sudPxVect, vector<Energy> & sudPyVect );

  /**
   * Along with the showering evolution --- going forward for
   * time-like (forward) evolution, and going backward for space-like
   * (backward) evolution --- the kinematical variables of the
   * branching products are calculated and updated from the knowledge
   * of the parent kinematics.  This method is used by the
   * ForwardShowerEvolver.  
   * ***ACHTUNG*** Might be extended to update colour connections as well.
   */
  virtual void updateChildren( const tShowerParticlePtr theParent, 
			       const ParticleVector theChildren );

  /**
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the 
   * KinematicsReconstructor.
   */
  virtual void updateParent( const tShowerParticlePtr theParent, 
			     const ParticleVector theChildren );

  /**
   * Update the kinematical data of a particle when a reconstruction
   * fixpoint was found. This will highly depend on the kind of
   * kinematics chosen and will be defined in the inherited concrete
   * classes. This method will be used by the KinematicsReconstructor.
   */
  virtual void updateLast( const tShowerParticlePtr theLast );

  /**
   * The method returns the mass of jet. 
   * It should be used only if isTheJetStartingPoint() 
   * is true, and only after the jet kinematics reconstruction.
   * (performed by the KinematicsReconstructor class).
   */
  virtual Energy jetMass();

private:

  /**
   * Private and non-existent assignment operator.
   */
  FS_QtildaShowerKinematics1to2 & operator=(const FS_QtildaShowerKinematics1to2 &);

};

}

#include "FS_QtildaShowerKinematics1to2.icc"

#endif /* HERWIG_FS_QtildaShowerKinematics1to2_H */
