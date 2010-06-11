// -*- C++ -*-
#ifndef HERWIG_SMHHHHVertex_H
#define HERWIG_SMHHHHVertex_H
//
// This is the declaration of the SMHHHHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/SSSSVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SMHHHHVertex class.
 *
 * @see \ref SMHHHHVertexInterfaces "The interfaces"
 * defined for SMHHHHVertex.
 */
  class SMHHHHVertex: public Helicity::SSSSVertex {

public:

  /**
   * The default constructor.
   */
  SMHHHHVertex();

  /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param part4 The ParticleData pointer for the fourth particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
			   tcPDPtr part3,tcPDPtr part4);

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SMHHHHVertex> initSMHHHHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHHHHVertex & operator=(const SMHHHHVertex &);

private:

  /**
   *  ratio of masses
   */
    double ratio_;

  /**
   *  The last value of the electroweak coupling calculated.
   */
  Complex couplast_;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 q2last_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMHHHHVertex. */
template <>
struct BaseClassTrait<Herwig::SMHHHHVertex,1> {
  /** Typedef of the first base class of SMHHHHVertex. */
  typedef Helicity::SSSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHHHHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMHHHHVertex>
  : public ClassTraitsBase<Herwig::SMHHHHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMHHHHVertex"; }
};

/** @endcond */

}

#endif /* HERWIG_SMHHHHVertex_H */
