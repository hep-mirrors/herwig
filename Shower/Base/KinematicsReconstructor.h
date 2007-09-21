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
#include "Herwig++/Shower/Nason/NasonTree.h"
#include "Evolver.fh"
#include "KinematicsReconstructor.fh"
#include <cassert>

namespace Herwig {

using namespace ThePEG;

  /**\ingroup Shower
   * Exception class
   * used to communicate failure of kinematics
   * reconstruction.
   */
  struct KinematicsReconstructionVeto {};


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
   * Given the ShowerTree for the shower from a hard process
   * the method does the reconstruction of the jets,
   * including the appropriate boosts (kinematics reshufflings)  
   * needed to conserve the total energy-momentum of the collision
   * and preserving the invariant mass and the rapidity of the 
   * hard subprocess system.
   */
  virtual bool reconstructHardJets(ShowerTreePtr hard,
				   map<tShowerProgenitorPtr,
				   pair<Energy,double> > pt) const=0;

  /**
   * Given the ShowerTree for a decay shower
   * the method does the reconstruction of the jets,
   * including the appropriate boosts (kinematics reshufflings)  
   * needed to conserve the total energy-momentum of the collision
   * and preserving the invariant mass and the rapidity of the 
   * hard subprocess system.
   */
  virtual bool reconstructDecayJets(ShowerTreePtr decay) const=0;
  //@}

  /**
   *  Methods to reconstruct the variables used to generate the
   *  shower given the particles produced.
   *  This is needed for the CKKW and Nason approaches
   */
  //@{
  /**
   *  Given the particles, with a history which we wish to interpret
   *  as a shower reconstruct the variables used to generate the 
   * shower for a decay process
   */
  virtual bool reconstructDecayShower(NasonTreePtr decay,EvolverPtr) const=0;

  /**
   *  Given the particles, with a history which we wish to interpret
   *  as a shower reconstruct the variables used to generate the shower
   *  for a hard process
   */
  virtual bool reconstructHardShower(NasonTreePtr hard,EvolverPtr) const=0;
  //@}

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
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static AbstractNoPIOClassDescription<KinematicsReconstructor> initKinematicsReconstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KinematicsReconstructor & operator=(const KinematicsReconstructor &);
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
  static string className() { return "Herwig::KinematicsReconstructor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * KinematicsReconstructor is implemented. It may also include several, space-separated,
   * libraries if the class KinematicsReconstructor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPI.so HwMPIPDF.so HwRemDecayer.so HwShower.so"; }
};

/** @endcond */

}

#include "KinematicsReconstructor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "KinematicsReconstructor.tcc"
#endif

#endif /* HERWIG_KinematicsReconstructor_H */
