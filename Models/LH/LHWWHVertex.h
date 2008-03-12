// -*- C++ -*-
#ifndef HERWIG_LHWWHVertex_H
#define HERWIG_LHWWHVertex_H
//
// This is the declaration of the LHWWHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "LHModel.h"
#include "LHWWHVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The LHWWHVertex class implements the couplings of two electroweak
 * gauge bosons to a Higgs boson in the Little Higgs model including the additional
 * heavy photon, Z and W bosons in the model and the triplet Higgs bosons.
 *
 * @see \ref LHWWHVertexInterfaces "The interfaces"
 * defined for LHWWHVertex.
 */
class LHWWHVertex: public Helicity::VVSVertex {

public:

  /**
   * The default constructor.
   */
  inline LHWWHVertex();
  
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
  static ClassDescription<LHWWHVertex> initLHWWHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHWWHVertex & operator=(const LHWWHVertex &);

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The last value of the electroweak coupling calculated.
   */
  Complex _couplast;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 _q2last;

  /**
   *  Couplings for the different interactions
   */
  vector<Energy> _coup;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LHWWHVertex. */
template <>
struct BaseClassTrait<Herwig::LHWWHVertex,1> {
  /** Typedef of the first base class of LHWWHVertex. */
  typedef Helicity::VVSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LHWWHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LHWWHVertex>
  : public ClassTraitsBase<Herwig::LHWWHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LHWWHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LHWWHVertex is implemented. It may also include several, space-separated,
   * libraries if the class LHWWHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLHModel.so"; }
};

/** @endcond */

}

#include "LHWWHVertex.icc"

#endif /* HERWIG_LHWWHVertex_H */
