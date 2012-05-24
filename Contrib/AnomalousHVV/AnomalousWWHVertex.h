// -*- C++ -*-
#ifndef HERWIG_AnomalousWWHVertex_H
#define HERWIG_AnomalousWWHVertex_H
//
// This is the declaration of the AnomalousWWHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/GeneralVVSVertex.h"

namespace Herwig {
using namespace ThePEG;

/**
 * The AnomalousWWHVertex class implements the option of an anomalous
 * coupling in the \f$h^0W^+W^+\f$ and \f$h^0Z^0Z^0\f$ vertices.
 * 
 * The options of using either the exact Standard Model, a CP-odd
 * \f[\frac1{\Lambda}J_1^\mu J_2^\nu\epsilon^{\mu\nu\alpha\beta}q_{1\alpha}q_{2\beta},\f]
 * or CP-even
 * \f[\frac1{\Lambda}J_1^\mu J_2^\nu
 *    \left[g_{\mu\nu}q_1 \cdot q_2-q_{1\nu}q_{2\mu}\right],\f]
 * form, where \f$\Lambda\f$ is the scale of the new operator.
 *
 * @see \ref AnomalousWWHVertexInterfaces "The interfaces"
 * defined for AnomalousWWHVertex.
 */
class AnomalousWWHVertex: public Helicity::GeneralVVSVertex {

public:

  /**
   * The default constructor.
   */
  AnomalousWWHVertex();

  /**
   * Calculate coupling.
   * @param q2 Scale at which to evaluate couplings
   * @param part1 ParticleDataPointer to first particle 
   * @param part2 ParticleDataPointer to second particle
   * @param part3 ParticleDataPointer to third particle 
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1, tcPDPtr part2,
                           tcPDPtr part3);

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
  static ClassDescription<AnomalousWWHVertex> initAnomalousWWHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AnomalousWWHVertex & operator=(const AnomalousWWHVertex &);

  /**
   *  Switch for the type of interaction
   */
  unsigned int interactionType_;

  /**
   *   The scale \f$\Lambda\f$
   */
  Energy Lambda_;

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The last value of the electroweak coupling calculated.
   */
  Complex couplast_;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 q2last_;

  /**
   *  The mass of the \f$W\f$ boson.
   */
  Energy mw_;

  /**
   *  The factor for the \f$Z\f$ vertex.
   */
  double zfact_;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AnomalousWWHVertex. */
template <>
struct BaseClassTrait<Herwig::AnomalousWWHVertex,1> {
  /** Typedef of the first base class of AnomalousWWHVertex. */
  typedef Helicity::GeneralVVSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AnomalousWWHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::AnomalousWWHVertex>
  : public ClassTraitsBase<Herwig::AnomalousWWHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::AnomalousWWHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * AnomalousWWHVertex is implemented. It may also include several, space-separated,
   * libraries if the class AnomalousWWHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "AnomalousHVV.so"; }
};

/** @endcond */

}

#endif /* HERWIG_AnomalousWWHVertex_H */
