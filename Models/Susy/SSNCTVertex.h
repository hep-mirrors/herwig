// -*- C++ -*-
#ifndef HERWIG_SSNCTVertex_H
#define HERWIG_SSNCTVertex_H
//
// This is the declaration of the SSNCTVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "MSSM.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SSNCTVertex class.
 *
 * @see \ref SSNCTVertexInterfaces "The interfaces"
 * defined for SSNCTVertex.
 */
class SSNCTVertex: public Helicity::FFSVertex {

public:

  /**
   * The default constructor.
   */
  SSNCTVertex();

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

  /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr part1,
                           tcPDPtr part2, tcPDPtr part3);

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
   * Initialize this object after the setup phase before saving and
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
  static ClassDescription<SSNCTVertex> initSSNCTVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSNCTVertex & operator=(const SSNCTVertex &);

private:

  /**
   *  High scale for the loops
   */
  Energy MX_;
  
  /**
   * Pointer to the neutralino mixing matrix
   */
  tMixingMatrixPtr nmix_;

  /**
   * \f$\sin(\theta_w)\f$
   */
  double sw_;

  /**
   * \f$\cos(\theta_w)\f$
   */
  double cw_;
  
  /**
   * Mass of the W
   */
  Energy mw_;

  /**
   * \f$\sin(\beta)\f$
   */
  double sb_;
  
  /**
   * \f$\cos(\beta)\f$
   */
  double cb_;

  /**
   * The scale at which the coupling was last evaluated. 
   */
  Energy2 q2last_;

  /**
   * The value of the normalisation when it was evaluated at _q2last 
   */
  Complex couplast_;
  
  /**
   * Store the value of the left coupling when it was last evaluated
   */
  Complex leftlast_;

  /**
   * Store the value of the right coupling when it was last evaluated
   */
  Complex rightlast_;

  /**
   * Store the id of the last neutralino to be evaluate
   */
  long idlast_;

  /**
   *  Mixing parameter
   */
  Complex epsilon_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSNCTVertex. */
template <>
struct BaseClassTrait<Herwig::SSNCTVertex,1> {
  /** Typedef of the first base class of SSNCTVertex. */
  typedef Helicity::FFSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSNCTVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSNCTVertex>
  : public ClassTraitsBase<Herwig::SSNCTVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SSNCTVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSNCTVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSNCTVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SSNCTVertex_H */
