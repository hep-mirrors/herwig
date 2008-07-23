// -*- C++ -*-
#ifndef HERWIG_LHTPFFPVertex_H
#define HERWIG_LHTPFFPVertex_H
//
// This is the declaration of the LHTPFFPVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "LHTPFFPVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the LHTPFFPVertex class.
 *
 * @see \ref LHTPFFPVertexInterfaces "The interfaces"
 * defined for LHTPFFPVertex.
 */
class LHTPFFPVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  LHTPFFPVertex();
  
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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<LHTPFFPVertex> initLHTPFFPVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHTPFFPVertex & operator=(const LHTPFFPVertex &);

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The charge of the Standard Model fermions.
   */
  vector<double> _charge;

  /**
   *  The last value of the coupling calculated.
   */
  Complex _couplast;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 _q2last;
  //@}

  /**
   *  Couplings of the fermion and T-fermions to the \f$A_H\f$
   */
  //@{
  /**
   *  Coupling of \f$dd_-A_H\f$
   */
  double _coupd;

  /**
   *  Coupling of \f$uu_-A_H\f$
   */
  double _coupu;

  /**
   *  Coupling of \f$ee_-A_H\f$
   */
  double _coupe;

  /**
   *  Coupling of \f$\nu\nu_-A_H\f$
   */
  double _coupnu;

  /**
   *  Coupling of \f$T_-T_+A_H\f$
   */
  double _tmtpL,_tmtpR;

  /**
   *  Coupling of \f$T_-tA_H\f$
   */
  double _tmtL,_tmtR;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LHTPFFPVertex. */
template <>
struct BaseClassTrait<Herwig::LHTPFFPVertex,1> {
  /** Typedef of the first base class of LHTPFFPVertex. */
  typedef Helicity::FFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LHTPFFPVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LHTPFFPVertex>
  : public ClassTraitsBase<Herwig::LHTPFFPVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LHTPFFPVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LHTPFFPVertex is implemented. It may also include several, space-separated,
   * libraries if the class LHTPFFPVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLHTPModel.so"; }
};

/** @endcond */

}

#endif /* HERWIG_LHTPFFPVertex_H */
