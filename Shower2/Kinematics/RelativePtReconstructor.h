// -*- C++ -*-
#ifndef HERWIG_RelativePtReconstructor_H
#define HERWIG_RelativePtReconstructor_H
//
// This is the declaration of the RelativePtReconstructor class.
//

#include "KinematicsReconstructor.h"
#include "RelativePtReconstructor.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the RelativePtReconstructor class.
 *
 * @see \ref RelativePtReconstructorInterfaces "The interfaces"
 * defined for RelativePtReconstructor.
 */
class RelativePtReconstructor: public KinematicsReconstructor {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline RelativePtReconstructor();

  /**
   * The copy constructor.
   */
  inline RelativePtReconstructor(const RelativePtReconstructor &);

  /**
   * The destructor.
   */
  virtual ~RelativePtReconstructor();
  //@}

public:

  /**
   * Given in input a vector of the particles which initiated the showers
   * the method does the reconstruction of such jets,
   * including the appropriate boosts (kinematics reshufflings)  
   * needed to conserve the total energy-momentum of the collision
   * and preserving the invariant mass and the rapidity of the 
   * hard subprocess system.
   */
  virtual bool reconstructHardJets(ShowerTreePtr hard);

protected:

  /**
   * Given the particle (ShowerParticle object) that 
   * originates a forward (time-like) jet, this method reconstructs the kinematics 
   * of the jet. That is, by starting from the final grand-children (which 
   * originates directly or indirectly from particleJetParent, 
   * and which don't have children), and moving "backwards" (in a physical
   * time picture), towards the particleJetParent, the 
   * ShowerKinematics objects associated with the various particles, 
   * which have been created during the showering, are now completed. 
   * In particular, at the end, we get the mass of the jet, which is the 
   * main information we want.
   * This methods returns false if there was no radiation or rescaling required
   */
  virtual bool reconstructTimeLikeJet(const tShowerParticlePtr particleJetParent,
				      unsigned int iopt);

  virtual void generateTimeLikeMomenta(const tShowerParticlePtr particleJetParent,
				       bool first);

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
   *  Calculate the off-shell mass of a particle
   */
  inline void calculateMass(const tShowerParticlePtr particleJetParent);

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
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<RelativePtReconstructor> initRelativePtReconstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RelativePtReconstructor & operator=(const RelativePtReconstructor &);

private:

  /**
   *  Option for the conversion between \f$\tilde{q}\f$ and mass.
   */
  unsigned int _massopt;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of RelativePtReconstructor. */
template <>
struct BaseClassTrait<Herwig::RelativePtReconstructor,1> {
  /** Typedef of the first base class of RelativePtReconstructor. */
  typedef Herwig::KinematicsReconstructor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the RelativePtReconstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::RelativePtReconstructor>
  : public ClassTraitsBase<Herwig::RelativePtReconstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::RelativePtReconstructor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * RelativePtReconstructor is implemented. It may also include several, space-separated,
   * libraries if the class RelativePtReconstructor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNewShower.so"; }
};

/** @endcond */

}

#include "RelativePtReconstructor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "RelativePtReconstructor.tcc"
#endif

#endif /* HERWIG_RelativePtReconstructor_H */
