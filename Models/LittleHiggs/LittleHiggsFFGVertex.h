// -*- C++ -*-
#ifndef HERWIG_LittleHiggsFFGVertex_H
#define HERWIG_LittleHiggsFFGVertex_H
//
// This is the declaration of the LittleHiggsFFGVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.fh"
#include "LittleHiggsFFGVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the LittleHiggsFFGVertex class.
 *
 * @see \ref LittleHiggsFFGVertexInterfaces "The interfaces"
 * defined for LittleHiggsFFGVertex.
 */
class LittleHiggsFFGVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  LittleHiggsFFGVertex();
  
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
  static NoPIOClassDescription<LittleHiggsFFGVertex> initLittleHiggsFFGVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LittleHiggsFFGVertex & operator=(const LittleHiggsFFGVertex &);
  
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
 *  base classes of LittleHiggsFFGVertex. */
template <>
struct BaseClassTrait<Herwig::LittleHiggsFFGVertex,1> {
  /** Typedef of the first base class of LittleHiggsFFGVertex. */
  typedef Helicity::FFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LittleHiggsFFGVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LittleHiggsFFGVertex>
  : public ClassTraitsBase<Herwig::LittleHiggsFFGVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LittleHiggsFFGVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LittleHiggsFFGVertex is implemented. It may also include several, space-separated,
   * libraries if the class LittleHiggsFFGVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLittleHiggsModel.so"; }
};

/** @endcond */

}

#include "LittleHiggsFFGVertex.icc"

#endif /* HERWIG_LittleHiggsFFGVertex_H */
