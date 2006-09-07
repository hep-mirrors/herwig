// -*- C++ -*-
#ifndef HERWIG_KinematicsReconstructor_H
#define HERWIG_KinematicsReconstructor_H
//
// This is the declaration of the KinematicsReconstructor class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ShowerParticle.h"
#include "ShowerProgenitor.h"
#include "ShowerTree.h"
#include "KinematicsReconstructor.fh"
#include <cassert>

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 * This class is responsible for the kinematical reconstruction 
 * after each showering step, and also for the necessary Lorentz boosts 
 * in order to preserve energy-momentum conservation in the overall collision,
 * and also the invariant mass and the rapidity of the hard subprocess system.
 * In the case of multi-step showering, there will be not unnecessary
 * kinematical reconstructions. 
 *
 * Notice:
 * - although we often use the term "jet" in either methods or variables names,
 *   or in comments, which could appear applicable only for QCD showering,
 *   there is indeed no "dynamics" represented in this class: only kinematics 
 *   is involved, as the name of this class remainds. Therefore it can be used
 *   for any kind of showers (QCD-,QED-,EWK-,... bremsstrahlung).
 * 
 * @see \ref KinematicsReconstructorInterfaces "The interfaces"
 * defined for KinematicsReconstructor.
 */
class KinematicsReconstructor: public Interfaced {

/**
 *  The Evolver is a friend to allow the  setting of _showerVariables
 */
friend class Evolver;

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
  virtual bool reconstructHardJets(ShowerTreePtr hard) const=0;

  /**
   * Given in input a vector of the particles which initiated the showers
   * the method does the reconstruction of such jets,
   * including the appropriate boosts (kinematics reshufflings)  
   * needed to conserve the total energy-momentum of the collision
   * and preserving the invariant mass and the rapidity of the 
   * hard subprocess system.
   */
  virtual bool reconstructDecayJets(ShowerTreePtr decay) const=0;
  //@}

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   *  Set/Get methods for the pointer to ShowerVariables
   */
  //@{
  
  /**
   *  Set the ShowerVariables
   */
  inline void showerVariables(ShowerVarsPtr);
  
  /**
   *  Get the ShowerVariables
   */
  inline ShowerVarsPtr showerVariables() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static AbstractClassDescription<KinematicsReconstructor> initKinematicsReconstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KinematicsReconstructor & operator=(const KinematicsReconstructor &);

private:

  /**
   *  Pointer to the ShowerVariables object
   */
  ShowerVarsPtr _showerVariables;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of KinematicsReconstructor. */
template <>
struct BaseClassTrait<Herwig::KinematicsReconstructor,1> {
  /** Typedef of the first base class of KinematicsReconstructor. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the KinematicsReconstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::KinematicsReconstructor>
  : public ClassTraitsBase<Herwig::KinematicsReconstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::KinematicsReconstructor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * KinematicsReconstructor is implemented. It may also include several, space-separated,
   * libraries if the class KinematicsReconstructor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "KinematicsReconstructor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "KinematicsReconstructor.tcc"
#endif

#endif /* HERWIG_KinematicsReconstructor_H */
