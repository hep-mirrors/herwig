// -*- C++ -*-
#ifndef HERWIG_QTildePTReconstructor_H
#define HERWIG_QTildePTReconstructor_H
//
// This is the declaration of the QTildePTReconstructor class.
//

#include "Herwig++/Shower/Base/KinematicsReconstructor.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup Shower
 *
 *
 * @see \ref QTildePTReconstructorInterfaces "The interfaces"
 * defined for QTildePTReconstructor.
 */
class QTildePTReconstructor: public KinematicsReconstructor {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline QTildePTReconstructor() : _massopt(0) {}
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
  virtual bool reconstructHardJets(ShowerTreePtr hard) const;

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
				      unsigned int iopt) const;

  /**
   *  Reconstruct the timelike momenta for a particle
   */
  virtual void generateTimeLikeMomenta(const tShowerParticlePtr particleJetParent,
				       bool first) const;

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
  inline void calculateMass(const tShowerParticlePtr particleJetParent) const
  {}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<QTildePTReconstructor> initQTildePTReconstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildePTReconstructor & operator=(const QTildePTReconstructor &);

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
 *  base classes of QTildePTReconstructor. */
template <>
struct BaseClassTrait<Herwig::QTildePTReconstructor,1> {
  /** Typedef of the first base class of QTildePTReconstructor. */
  typedef Herwig::KinematicsReconstructor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QTildePTReconstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QTildePTReconstructor>
  : public ClassTraitsBase<Herwig::QTildePTReconstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QTildePTReconstructor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QTildePTReconstructor is implemented. It may also include several, space-separated,
   * libraries if the class QTildePTReconstructor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPIPDF.so HwRemDecayer.so HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_QTildePTReconstructor_H */
