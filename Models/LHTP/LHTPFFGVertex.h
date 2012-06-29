// -*- C++ -*-
#ifndef HERWIG_LHTPFFGVertex_H
#define HERWIG_LHTPFFGVertex_H
//
// This is the declaration of the LHTPFFGVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "LHTPFFGVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The LHTPFFGVertex class implements the coupling of the 
 * gluon to the coloured fermions, the SM quarks, the extra top-like quark
 * and the T-parity odd quarks of the Little Higgs model with T-parity.
 *
 * @see \ref LHTPFFGVertexInterfaces "The interfaces"
 * defined for LHTPFFGVertex.
 */
class LHTPFFGVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  LHTPFFGVertex();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<LHTPFFGVertex> initLHTPFFGVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHTPFFGVertex & operator=(const LHTPFFGVertex &);

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The last value of the strong coupling calculated.
   */
  Complex _couplast;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 _q2last;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LHTPFFGVertex. */
template <>
struct BaseClassTrait<Herwig::LHTPFFGVertex,1> {
  /** Typedef of the first base class of LHTPFFGVertex. */
  typedef Helicity::FFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LHTPFFGVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LHTPFFGVertex>
  : public ClassTraitsBase<Herwig::LHTPFFGVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LHTPFFGVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LHTPFFGVertex is implemented. It may also include several, space-separated,
   * libraries if the class LHTPFFGVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLHTPModel.so"; }
};

/** @endcond */

}

#endif /* HERWIG_LHTPFFGVertex_H */
