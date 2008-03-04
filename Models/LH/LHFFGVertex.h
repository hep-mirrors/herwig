// -*- C++ -*-
#ifndef HERWIG_LHFFGVertex_H
#define HERWIG_LHFFGVertex_H
//
// This is the declaration of the LHFFGVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.fh"
#include "LHFFGVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The LHFFGVertex class implements the couplings of the 
 * quarks and additional heavy top to the gluon in the Little Higgs
 * model.
 *
 * @see \ref LHFFGVertexInterfaces "The interfaces"
 * defined for LHFFGVertex.
 */
class LHFFGVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  LHFFGVertex();
  
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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<LHFFGVertex> initLHFFGVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHFFGVertex & operator=(const LHFFGVertex &);
  
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
 *  base classes of LHFFGVertex. */
template <>
struct BaseClassTrait<Herwig::LHFFGVertex,1> {
  /** Typedef of the first base class of LHFFGVertex. */
  typedef Helicity::FFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LHFFGVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LHFFGVertex>
  : public ClassTraitsBase<Herwig::LHFFGVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LHFFGVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LHFFGVertex is implemented. It may also include several, space-separated,
   * libraries if the class LHFFGVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLHModel.so"; }
};

/** @endcond */

}

#include "LHFFGVertex.icc"

#endif /* HERWIG_LHFFGVertex_H */
