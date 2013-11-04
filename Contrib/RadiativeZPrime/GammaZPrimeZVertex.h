// -*- C++ -*-
#ifndef RADIATIVEZPRIME_GammaZPrimeZVertex_H
#define RADIATIVEZPRIME_GammaZPrimeZVertex_H
//
// This is the declaration of the GammaZPrimeZVertex class.
//

#include "AnomalousVVVVertex.h"

namespace RadiativeZPrime {

using namespace ThePEG;

/**
 * Here is the documentation of the GammaZPrimeZVertex class.
 *
 * @see \ref GammaZPrimeZVertexInterfaces "The interfaces"
 * defined for GammaZPrimeZVertex.
 */
class GammaZPrimeZVertex: public AnomalousVVVVertex {

public:

  /**
   * The default constructor.
   */
  GammaZPrimeZVertex();

  /**
   * Calculate the couplings. This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

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
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<GammaZPrimeZVertex> initGammaZPrimeZVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GammaZPrimeZVertex & operator=(const GammaZPrimeZVertex &);

private:

  /**
   *   Storage of the coupling
   */
  double _coup;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GammaZPrimeZVertex. */
template <>
struct BaseClassTrait<RadiativeZPrime::GammaZPrimeZVertex,1> {
  /** Typedef of the first base class of GammaZPrimeZVertex. */
  typedef RadiativeZPrime::AnomalousVVVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GammaZPrimeZVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<RadiativeZPrime::GammaZPrimeZVertex>
  : public ClassTraitsBase<RadiativeZPrime::GammaZPrimeZVertex> {
  /** Return a platform-independent class name */
  static string className() { return "RadiativeZPrime::GammaZPrimeZVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * GammaZPrimeZVertex is implemented. It may also include several, space-separated,
   * libraries if the class GammaZPrimeZVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "RadiativeZPrime.so"; }
};

/** @endcond */

}

#endif /* RADIATIVEZPRIME_GammaZPrimeZVertex_H */
