// -*- C++ -*-
#ifndef HERWIG_LHWWWVertex_H
#define HERWIG_LHWWWVertex_H
//
// This is the declaration of the LHWWWVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"
#include "LHModel.h"
#include "LHWWWVertex.fh"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::Direction;

/**
 * The LHWWWVertex class implements the triple boson coupling in the Little
 * Higgs model for the electroweak bosons of the Standard Model and the
 * additional \f$A_H\f$, \f$Z_H\f$ and \f$W_H^\pm\f$ bosons.
 *
 * @see \ref LHWWWVertexInterfaces "The interfaces"
 * defined for LHWWWVertex.
 */
class LHWWWVertex: public Helicity::VVVVertex {

public:

  /**
   * The default constructor.
   */
  inline LHWWWVertex();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param d1 The direction for the first  particle.
   * @param d2 The direction for the second particle.
   * @param d3 The direction for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3,
			   Direction d1,Direction d2, Direction d3);

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
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<LHWWWVertex> initLHWWWVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHWWWVertex & operator=(const LHWWWVertex &);

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
 *  base classes of LHWWWVertex. */
template <>
struct BaseClassTrait<Herwig::LHWWWVertex,1> {
  /** Typedef of the first base class of LHWWWVertex. */
  typedef Helicity::VVVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LHWWWVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LHWWWVertex>
  : public ClassTraitsBase<Herwig::LHWWWVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LHWWWVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LHWWWVertex is implemented. It may also include several, space-separated,
   * libraries if the class LHWWWVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLHModel.so"; }
};

/** @endcond */

}

#include "LHWWWVertex.icc"

#endif /* HERWIG_LHWWWVertex_H */
