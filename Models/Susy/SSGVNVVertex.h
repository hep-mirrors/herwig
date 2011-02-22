// -*- C++ -*-
#ifndef HERWIG_SSGVNVVertex_H
#define HERWIG_SSGVNVVertex_H
//
// This is the declaration of the SSGVNVVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/RFVVertex.h"
#include "MixingMatrix.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The SSGVNVVertex class implements the coupling of the gravitino
 * to a gaugino and the assoicated gauge boson.
 *
 * @see \ref SSGVNVVertexInterfaces "The interfaces"
 * defined for SSGVNVVertex.
 */
class SSGVNVVertex: public Helicity::RFVVertex {

public:

  /**
   * The default constructor.
   */
  SSGVNVVertex();

  /**
   * Calculate the couplings. This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
			   tcPDPtr part2,tcPDPtr part3);

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
  static ClassDescription<SSGVNVVertex> initSSGVNVVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSGVNVVertex & operator=(const SSGVNVVertex &);

private:

  /**
   *  \f$\sin\theta_W\f$
   */
  double sw_;

  /**
   *  \f$\sin\theta_W\f$
   */
  double cw_;

  /**
   * The value of \f$\sin\beta\f$ 
   */
  double sb_;

  /**
   * The value of \f$\cos\beta\f$ 
   */
  double cb_;

  /**
   *  The Z mass
   */  
  Energy mz_;
  
  /**
   * Pointer to the neutralino mixing matrix
   */
  tMixingMatrixPtr nmix_;

  /**
   *  The Planck mass
   */
  Energy MPlanck_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSGVNVVertex. */
template <>
struct BaseClassTrait<Herwig::SSGVNVVertex,1> {
  /** Typedef of the first base class of SSGVNVVertex. */
  typedef Helicity::RFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSGVNVVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSGVNVVertex>
  : public ClassTraitsBase<Herwig::SSGVNVVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SSGVNVVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSGVNVVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSGVNVVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "SSGVNVVertex.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SSGVNVVertex_H */
