// -*- C++ -*-
#ifndef HERWIG_LHWHHVertex_H
#define HERWIG_LHWHHVertex_H
//
// This is the declaration of the LHWHHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the LHWHHVertex class.
 *
 * @see \ref LHWHHVertexInterfaces "The interfaces"
 * defined for LHWHHVertex.
 */
class LHWHHVertex: public Helicity::VSSVertex {

public:

  /**
   * The default constructor.
   */
  LHWHHVertex();

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
   * Calculate the coupling for the vertex
   * @param q2 The scale to at which evaluate the coupling.
   * @param particle1 The first particle in the vertex.
   * @param particle2 The second particle in the vertex.
   * @param particle3 The third particle in the vertex.
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			   tcPDPtr particle3);

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
  static ClassDescription<LHWHHVertex> initLHWHHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHWHHVertex & operator=(const LHWHHVertex &);

private:

  /**
   * The value of the coupling when last evaluated
   */
  Complex couplast_;
  
  /**
   * The scale at which the coupling  was last evaluated.
   */
  Energy2 q2last_;

  /**
   *  Couplings 
   */
  vector<Complex> coup_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LHWHHVertex. */
template <>
struct BaseClassTrait<Herwig::LHWHHVertex,1> {
  /** Typedef of the first base class of LHWHHVertex. */
  typedef Helicity::VSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LHWHHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LHWHHVertex>
  : public ClassTraitsBase<Herwig::LHWHHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LHWHHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LHWHHVertex is implemented. It may also include several, space-separated,
   * libraries if the class LHWHHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLHModel.so"; }
};

/** @endcond */

}

#endif /* HERWIG_LHWHHVertex_H */
