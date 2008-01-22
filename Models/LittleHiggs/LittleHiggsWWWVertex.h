// -*- C++ -*-
#ifndef HERWIG_LittleHiggsWWWVertex_H
#define HERWIG_LittleHiggsWWWVertex_H
//
// This is the declaration of the LittleHiggsWWWVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"
#include "LittleHiggsModel.h"
#include "LittleHiggsWWWVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the LittleHiggsWWWVertex class.
 *
 * @see \ref LittleHiggsWWWVertexInterfaces "The interfaces"
 * defined for LittleHiggsWWWVertex.
 */
class LittleHiggsWWWVertex: public Helicity::VVVVertex {

public:

  /**
   * The default constructor.
   */
  inline LittleHiggsWWWVertex();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

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
  static ClassDescription<LittleHiggsWWWVertex> initLittleHiggsWWWVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LittleHiggsWWWVertex & operator=(const LittleHiggsWWWVertex &);

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The correction factors for the different interacting particles
   */
  vector<double> _corr;

  /**
   *  The last value of the electroweak coupling calculated.
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
 *  base classes of LittleHiggsWWWVertex. */
template <>
struct BaseClassTrait<Herwig::LittleHiggsWWWVertex,1> {
  /** Typedef of the first base class of LittleHiggsWWWVertex. */
  typedef Helicity::VVVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LittleHiggsWWWVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LittleHiggsWWWVertex>
  : public ClassTraitsBase<Herwig::LittleHiggsWWWVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LittleHiggsWWWVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LittleHiggsWWWVertex is implemented. It may also include several, space-separated,
   * libraries if the class LittleHiggsWWWVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLittleHiggsModel.so"; }
};

/** @endcond */

}

#include "LittleHiggsWWWVertex.icc"

#endif /* HERWIG_LittleHiggsWWWVertex_H */
