// -*- C++ -*-
#ifndef HERWIG_IS_QtildaShowerKinematics1to2_H
#define HERWIG_IS_QtildaShowerKinematics1to2_H
//
// This is the declaration of the IS_QtildaShowerKinematics1to2 class.

#include "QtildaShowerKinematics1to2.h"


namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  
 *  This (concrete) class provides the specific Intial State shower
 *  kinematics information.
 *
 *  @see QtildaShowerKinematics1to2
 *  @see FS_QtildaShowerKinematics1to2
 *  @see KinematicsReconstructor
 */
class IS_QtildaShowerKinematics1to2: public QtildaShowerKinematics1to2 {

public:

  /**
   * Standard ctors and dtor.
   */
  inline IS_QtildaShowerKinematics1to2();
  inline IS_QtildaShowerKinematics1to2(const IS_QtildaShowerKinematics1to2 &);
  inline IS_QtildaShowerKinematics1to2( const Lorentz5Momentum & p, 
					const Lorentz5Momentum & n );
  virtual ~IS_QtildaShowerKinematics1to2();

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
   * Update the parent Kinematics from the knowledge of the kinematics
   * of the children. This method will be used by the KinematicsReconstructor.
   */
  virtual void updateParent( tCollecShoKinPtr & shoKinChildren );

  virtual void updateChildren( const tShowerParticlePtr theParent, 
			       const ShowerParticleVector theChildren );
  virtual void updateParent( const tShowerParticlePtr theParent, 
			     const ParticleVector theChildren );
  virtual void updateLast( const tShowerParticlePtr theLast );

  /**
   * The method returns the mass of jet. 
   * It should be used only if isTheJetStartingPoint() is true, 
   * and only after the jet kinematics reconstruction.
   * (performed by the KinematicsReconstructor class).
   */
  virtual Energy jetMass();

private:

  /**
   * Private and non-existent assignment operator.
   */
  IS_QtildaShowerKinematics1to2 & operator=(const IS_QtildaShowerKinematics1to2 &);

};

}

#include "IS_QtildaShowerKinematics1to2.icc"

#endif /* HERWIG_IS_QtildaShowerKinematics1to2_H */
