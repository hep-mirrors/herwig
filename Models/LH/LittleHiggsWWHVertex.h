// -*- C++ -*-
#ifndef HERWIG_LittleHiggsWWHVertex_H
#define HERWIG_LittleHiggsWWHVertex_H
//
// This is the declaration of the LittleHiggsWWHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "LittleHiggsModel.h"
#include "LittleHiggsWWHVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The LittleHiggsWWHVertex class implements the couplings of two electroweak
 * gauge bosons to a Higgs boson in the Little Higgs model including the additional
 * heavy photon, Z and W bosons in the model and the triplet Higgs bosons.
 *
 * @see \ref LittleHiggsWWHVertexInterfaces "The interfaces"
 * defined for LittleHiggsWWHVertex.
 */
class LittleHiggsWWHVertex: public Helicity::VVSVertex {

public:

  /**
   * The default constructor.
   */
  inline LittleHiggsWWHVertex();
  
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
  static ClassDescription<LittleHiggsWWHVertex> initLittleHiggsWWHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LittleHiggsWWHVertex & operator=(const LittleHiggsWWHVertex &);

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
 *  base classes of LittleHiggsWWHVertex. */
template <>
struct BaseClassTrait<Herwig::LittleHiggsWWHVertex,1> {
  /** Typedef of the first base class of LittleHiggsWWHVertex. */
  typedef Helicity::VVSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LittleHiggsWWHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LittleHiggsWWHVertex>
  : public ClassTraitsBase<Herwig::LittleHiggsWWHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LittleHiggsWWHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LittleHiggsWWHVertex is implemented. It may also include several, space-separated,
   * libraries if the class LittleHiggsWWHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLHModel.so"; }
};

/** @endcond */

}

#include "LittleHiggsWWHVertex.icc"

#endif /* HERWIG_LittleHiggsWWHVertex_H */
