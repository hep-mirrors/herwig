// -*- C++ -*-
#ifndef RADIATIVEZPRIME_FFZPrimeVertex_H
#define RADIATIVEZPRIME_FFZPrimeVertex_H
//
// This is the declaration of the FFZPrimeVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"

namespace RadiativeZPrime {

using namespace ThePEG;

/**
 * Here is the documentation of the FFZPrimeVertex class.
 *
 * @see \ref FFZPrimeVertexInterfaces "The interfaces"
 * defined for FFZPrimeVertex.
 */
class FFZPrimeVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  inline FFZPrimeVertex();
  
  /**
   * Calculate the couplings. 
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
  
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<FFZPrimeVertex> initFFZPrimeVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FFZPrimeVertex & operator=(const FFZPrimeVertex &);


  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The left couplings of the Standard Model fermions.
   */
  vector<double> _gl;

  /**
   *  The right couplings of the Standard Model fermions.
   */
  vector<double> _gr;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FFZPrimeVertex. */
template <>
struct BaseClassTrait<RadiativeZPrime::FFZPrimeVertex,1> {
  /** Typedef of the first base class of FFZPrimeVertex. */
  typedef Helicity::FFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FFZPrimeVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<RadiativeZPrime::FFZPrimeVertex>
  : public ClassTraitsBase<RadiativeZPrime::FFZPrimeVertex> {
  /** Return a platform-independent class name */
  static string className() { return "RadiativeZPrime::FFZPrimeVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FFZPrimeVertex is implemented. It may also include several, space-separated,
   * libraries if the class FFZPrimeVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "RadiativeZPrime.so"; }
};

/** @endcond */

}

#endif /* RADIATIVEZPRIME_FFZPrimeVertex_H */
