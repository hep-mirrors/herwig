// -*- C++ -*-
#ifndef HERWIG_GravitonMassGenerator_H
#define HERWIG_GravitonMassGenerator_H
//
// This is the declaration of the GravitonMassGenerator class.
//

#include "Herwig/PDT/GenericMassGenerator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the GravitonMassGenerator class.
 *
 * @see \ref GravitonMassGeneratorInterfaces "The interfaces"
 * defined for GravitonMassGenerator.
 */
class GravitonMassGenerator: public GenericMassGenerator {

public:

  /**
   * The default constructor.
   */
  GravitonMassGenerator();

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

protected:

  /**
   * Return a mass with the weight using the specified limits.
   * @param part The particle data pointer of the particle.
   * @param low The lower limit on the particle's mass.
   * @param upp The upper limit on the particle's mass.
   * @param wgt The weight for this mass.
   * @param shape The type of shape to use
   * @param r   The random number used for the weight
   * @return The mass of the particle instance.
   */
  virtual Energy mass(double & wgt, const ParticleData & part,
		      const Energy low,const Energy upp, int shape,
		      double r=UseRandom::rnd()) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GravitonMassGenerator & operator=(const GravitonMassGenerator &) = delete;

private:

  /**
   *  prefactor
   */
  double prefactor_;

  /**
   *  Number of extra dimensions
   */
  unsigned int delta_;

  /**
   *  d-dimensional Planck mass
   */
  Energy md_;

  /**
   *  Minimum mass cut to avoid numerical problems
   */
  Energy mMin_;

};

}

#endif /* HERWIG_GravitonMassGenerator_H */
