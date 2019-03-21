// -*- C++ -*-
#ifndef HERWIG_SSGVNHVertex_H
#define HERWIG_SSGVNHVertex_H
//
// This is the declaration of the SSGVNHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/RFSVertex.h"
#include "MixingMatrix.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The SSGVNHVertex class implements the coupling of the gravitino to
 * a neutralino and Higgs boson.
 *
 * @see \ref SSGVNHVertexInterfaces "The interfaces"
 * defined for SSGVNHVertex.
 */
class SSGVNHVertex: public Helicity::RFSVertex {

public:

  /**
   * The default constructor.
   */
  SSGVNHVertex();

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
  static ClassDescription<SSGVNHVertex> initSSGVNHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSGVNHVertex & operator=(const SSGVNHVertex &) = delete;

private:

  /**
   * The value of \f$\sin\alpha\f$ 
   */
  double sa_;

  /**
   * The value of \f$\sin\beta\f$ 
   */
  double sb_;

  /**
   * The value of \f$\cos\alpha\f$ 
   */
  double ca_;

  /**
   * The value of \f$\cos\beta\f$ 
   */
  double cb_;  
  
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
 *  base classes of SSGVNHVertex. */
template <>
struct BaseClassTrait<Herwig::SSGVNHVertex,1> {
  /** Typedef of the first base class of SSGVNHVertex. */
  typedef Helicity::RFSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSGVNHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSGVNHVertex>
  : public ClassTraitsBase<Herwig::SSGVNHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SSGVNHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSGVNHVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSGVNHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SSGVNHVertex_H */
