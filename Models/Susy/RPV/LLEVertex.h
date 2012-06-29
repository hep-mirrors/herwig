// -*- C++ -*-
#ifndef HERWIG_LLEVertex_H
#define HERWIG_LLEVertex_H
//
// This is the declaration of the LLEVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "RPV.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the LLEVertex class.
 *
 * @see \ref LLEVertexInterfaces "The interfaces"
 * defined for LLEVertex.
 */
class LLEVertex: public FFSVertex {

public:

  /**
   * The default constructor.
   */
  LLEVertex();

  /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr part1,
                           tcPDPtr part2, tcPDPtr part3);

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<LLEVertex> initLLEVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LLEVertex & operator=(const LLEVertex &);

private:

  /**
   *  Coupling
   */
  vector<vector<vector<double > > > lambda_;

  /**
   * Pointer to the stau mixing matrix
   */
  tMixingMatrixPtr stau_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LLEVertex. */
template <>
struct BaseClassTrait<Herwig::LLEVertex,1> {
  /** Typedef of the first base class of LLEVertex. */
  typedef FFSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LLEVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LLEVertex>
  : public ClassTraitsBase<Herwig::LLEVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LLEVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LLEVertex is implemented. It may also include several, space-separated,
   * libraries if the class LLEVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "LLEVertex.so"; }
};

/** @endcond */

}

#endif /* HERWIG_LLEVertex_H */
