// -*- C++ -*-
#ifndef HERWIG_FortranReconstructor_H
#define HERWIG_FortranReconstructor_H
//
// This is the declaration of the FortranReconstructor class.
//

#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "Herwig++/Shower/ShowerProgenitor.h"
#include "Herwig++/Shower/ShowerTree.h"
#include "FortranReconstructor.fh"
#include <cassert>

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 * This class is responsible for the kinematical reconstruction 
 * after each showering step for the FORTRAN HERWIG shower.
 * 
 * @see ShowerParticle
 * @see ShowerKinematics
 * @see \ref FortranReconstructorInterfaces "The interfaces"
 * defined for FortranReconstructor.
 */
class FortranReconstructor: public KinematicsReconstructor {

public:

  /**
   *  Methods to reconstruct the kinematics of a scattering or decay process
   */
  //@{
  /**
   * Given in input a vector of the particles which initiated the showers
   * the method does the reconstruction of such jets,
   * including the appropriate boosts (kinematics reshufflings)  
   * needed to conserve the total energy-momentum of the collision
   * and preserving the invariant mass and the rapidity of the 
   * hard subprocess system.
   */
  virtual bool reconstructHardJets(ShowerTreePtr hard) const;

  /**
   * Given in input a vector of the particles which initiated the showers
   * the method does the reconstruction of such jets,
   * including the appropriate boosts (kinematics reshufflings)  
   * needed to conserve the total energy-momentum of the collision
   * and preserving the invariant mass and the rapidity of the 
   * hard subprocess system.
   */
  virtual bool reconstructDecayJets(ShowerTreePtr decay) const;
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<FortranReconstructor> initFortranReconstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FortranReconstructor & operator=(const FortranReconstructor &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FortranReconstructor. */
template <>
struct BaseClassTrait<Herwig::FortranReconstructor,1> {
  /** Typedef of the first base class of FortranReconstructor. */
  typedef Herwig::KinematicsReconstructor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FortranReconstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FortranReconstructor>
  : public ClassTraitsBase<Herwig::FortranReconstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::FortranReconstructor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FortranReconstructor is implemented. It may also include several, space-separated,
   * libraries if the class FortranReconstructor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "FortranReconstructor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FortranReconstructor.tcc"
#endif

#endif /* HERWIG_FortranReconstructor_H */
