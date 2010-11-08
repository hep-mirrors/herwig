// -*- C++ -*-
#ifndef HERWIG_UDDVertex_H
#define HERWIG_UDDVertex_H
//
// This is the declaration of the UDDVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "RPV.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the UDDVertex class.
 *
 * @see \ref UDDVertexInterfaces "The interfaces"
 * defined for UDDVertex.
 */
class UDDVertex: public FFSVertex {

public:

  /**
   * The default constructor.
   */
  UDDVertex();

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
  static ClassDescription<UDDVertex> initUDDVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  UDDVertex & operator=(const UDDVertex &);

private:

  /**
   *  Coupling
   */
  vector<vector<vector<double > > > lambda_;

  /**
   * Pointer to the stop mixing matrix
   */
  tMixingMatrixPtr stop_;

  /**
   * Pointer to the sbottom mixing matrix
   */
  tMixingMatrixPtr sbot_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of UDDVertex. */
template <>
struct BaseClassTrait<Herwig::UDDVertex,1> {
  /** Typedef of the first base class of UDDVertex. */
  typedef FFSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the UDDVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::UDDVertex>
  : public ClassTraitsBase<Herwig::UDDVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::UDDVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * UDDVertex is implemented. It may also include several, space-separated,
   * libraries if the class UDDVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "UDDVertex.so"; }
};

/** @endcond */

}

#endif /* HERWIG_UDDVertex_H */
