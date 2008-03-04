// -*- C++ -*-
#ifndef HERWIG_LHTPFFZVertex_H
#define HERWIG_LHTPFFZVertex_H
//
// This is the declaration of the LHTPFFZVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "LHTPFFZVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the LHTPFFZVertex class.
 *
 * @see \ref LHTPFFZVertexInterfaces "The interfaces"
 * defined for LHTPFFZVertex.
 */
class LHTPFFZVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  inline LHTPFFZVertex();
  
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
  static ClassDescription<LHTPFFZVertex> initLHTPFFZVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHTPFFZVertex & operator=(const LHTPFFZVertex &);

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The left couplings of the Standard Model fermions.
   */
  vector<double> _gl;

  /**
   *  The right couplings of the Standard Model fermions.
   */
  vector<double> _gr;

  /**
   *  The left couplings to the Z for the extended top sector
   */
  vector<double> _tl;

  /**
   *  The right couplings to the Z for the extended top sector
   */
  vector<double> _tr;

  /**
   *  Coupling of \f$dd_-Z_H\f$
   */
  double _coupd;

  /**
   *  Coupling of \f$uu_-Z_H\f$
   */
  double _coupu;

  /**
   *  Coupling of \f$ee_-Z_H\f$
   */
  double _coupe;

  /**
   *  Coupling of \f$\nu\nu_-Z_H\f$
   */
  double _coupnu;

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
 *  base classes of LHTPFFZVertex. */
template <>
struct BaseClassTrait<Herwig::LHTPFFZVertex,1> {
  /** Typedef of the first base class of LHTPFFZVertex. */
  typedef Helicity::FFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LHTPFFZVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LHTPFFZVertex>
  : public ClassTraitsBase<Herwig::LHTPFFZVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LHTPFFZVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LHTPFFZVertex is implemented. It may also include several, space-separated,
   * libraries if the class LHTPFFZVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLHTP.so"; }
};

/** @endcond */

}

#include "LHTPFFZVertex.icc"

#endif /* HERWIG_LHTPFFZVertex_H */
