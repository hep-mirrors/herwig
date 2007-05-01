// -*- C++ -*-
#ifndef HERWIG_SSGFSVertex_H
#define HERWIG_SSGFSVertex_H
//
// This is the declaration of the SSGFSVertex class.
//

#include "Herwig++/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Models/Susy/SusyBase.h"
#include "SSGFSVertex.fh"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/**
 * This is the implementation of the coupling of a gluino to a squark-
 * quark pair. It inherits from FFSVertex and implements the setCoupling
 * method.
 *
 * @see FFSVertex
 */
class SSGFSVertex: public FFSVertex {

public:

  /**
   * The default constructor.
   */
  SSGFSVertex();

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

  /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param ioff Integer giving the off-shell particle
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
                           tcPDPtr part2,tcPDPtr part3, int ioff);

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
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  //@}

  /**
   * Pointer to the SusyBase object
   */
  tSusyBasePtr _theSS;
  
  /**
   * Pointer to the stop mixing matrix
   */
  tMixingMatrixPtr _stop;

  /**
   * Pointer to the _sbottom mixing matrix
   */
  tMixingMatrixPtr _sbottom;
  
  /**
   * The scale at which the coupling was last evaluated;
   */
  Energy2 _q2last;

  /**
   * The value of the coupling when it was last evaluated
   */
  Complex _couplast;

  /**
   * The id of the sm fermion for which the coupling was evaluated
   */
  long _id1last;
  
  /**
   * The id of the scalar for which the coupling was evaluated
   */
  long _id2last;

  /**
   * The value of the left coupling when it was last evaluated
   */
  Complex _leftlast;

  /**
   * The value of the right coupling when it was last evaluated
   */
  Complex _rightlast;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SSGFSVertex> initSSGFSVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSGFSVertex & operator=(const SSGFSVertex &);

};
}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSGFSVertex. */
template <>
struct BaseClassTrait<Herwig::Helicity::SSGFSVertex,1> {
  /** Typedef of the first base class of SSGFSVertex. */
  typedef Herwig::Helicity::FFSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSGFSVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Helicity::SSGFSVertex>
  : public ClassTraitsBase<Herwig::Helicity::SSGFSVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SSGFSVertex"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SSGFSVertex class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwSusy.so HwSusyVertex.so"; }
};

/** @endcond */

}

#include "SSGFSVertex.icc"

#endif /* HERWIG_SSGFSVertex_H */
