// -*- C++ -*-
#ifndef HERWIG_SMZDecayerPOWHEG_H
#define HERWIG_SMZDecayerPOWHEG_H
//
// This is the declaration of the SMZDecayerPOWHEG class.
//

#include "SMZDecayer.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SMZDecayerPOWHEG class.
 *
 * @see \ref SMZDecayerPOWHEGInterfaces "The interfaces"
 * defined for SMZDecayerPOWHEG.
 */
class SMZDecayerPOWHEG: public SMZDecayer {

public:

  /**
   * The default constructor.
   */
  SMZDecayerPOWHEG() : contrib_(1), zPow_(0.5), yPow_(0.9)
  {}

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int ichan, const Particle & part,
		     const ParticleVector & decay,MEOption meopt) const;

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

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SMZDecayerPOWHEG> initSMZDecayerPOWHEG;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMZDecayerPOWHEG & operator=(const SMZDecayerPOWHEG &);

private:

  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Phase-space sampling for z
   */
  double zPow_;

  /**
   *  Phase-space sampling for y
   */
  double yPow_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMZDecayerPOWHEG. */
template <>
struct BaseClassTrait<Herwig::SMZDecayerPOWHEG,1> {
  /** Typedef of the first base class of SMZDecayerPOWHEG. */
  typedef SMZDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMZDecayerPOWHEG class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMZDecayerPOWHEG>
  : public ClassTraitsBase<Herwig::SMZDecayerPOWHEG> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMZDecayerPOWHEG"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SMZDecayerPOWHEG is implemented. It may also include several, space-separated,
   * libraries if the class SMZDecayerPOWHEG depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "SMZDecayerPOWHEG.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SMZDecayerPOWHEG_H */
