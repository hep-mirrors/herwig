// -*- C++ -*-
#ifndef HERWIG_MEvv2ss_H
#define HERWIG_MEvv2ss_H
//
// This is the declaration of the MEvv2ss class.
//

#include "GeneralHardME.h"
#include "MEvv2ss.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This is the implementation of the matrix element for the process
 * vector-vector to scalar-scalar. It inherits from GeneralHardME and
 * implements the required virtual functions.
 *
 * @see \ref MEff2ffInterfaces "The Interfaces"
 * defined for MEff2ff.
 * @see GeneralHardME
 */
class MEvv2ss: public GeneralHardME {

public:

  /**
   * The default constructor.
   */
  inline MEvv2ss();

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
  //@}

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEvv2ss> initMEvv2ss;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEvv2ss & operator=(const MEvv2ss &);
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEvv2ss. */
template <>
struct BaseClassTrait<Herwig::MEvv2ss,1> {
  /** Typedef of the first base class of MEvv2ss. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEvv2ss class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEvv2ss>
  : public ClassTraitsBase<Herwig::MEvv2ss> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEvv2ss"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEvv2ss is implemented. It may also include several, space-separated,
   * libraries if the class MEvv2ss depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libGeneralHardME.so"; }
};

/** @endcond */

}

#include "MEvv2ss.icc"

#endif /* HERWIG_MEvv2ss_H */
