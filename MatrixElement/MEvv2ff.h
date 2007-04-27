// -*- C++ -*-
#ifndef HERWIG_MEvv2ff_H
#define HERWIG_MEvv2ff_H
//
// This is the declaration of the MEvv2ff class.
//

#include "GeneralHardME.h"
#include "MEvv2ff.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class is designed to implement the matrix element for the 
 * \f$2 \rightarrow 2\f$ process vector-vector to fermion-antifermion pair. It
 * inherits from GeneralHardME and implements the me2() virtual function.
 *
 * @see \ref MEvv2ffInterfaces "The Interfaces"
 * defined for MEvv2ff.
 * @see GeneralHardME
 * 
 */
class MEvv2ff: public GeneralHardME {

public:

  /**
   * The default constructor.
   */
  inline MEvv2ff();

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
  static ClassDescription<MEvv2ff> initMEvv2ff;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEvv2ff & operator=(const MEvv2ff &);
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/// \if TRAITSPECIALIZATIONS

/** This template specialization informs ThePEG about the
 *  base classes of MEvv2ff. */
template <>
struct BaseClassTrait<Herwig::MEvv2ff,1> {
  /** Typedef of the first base class of MEvv2ff. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEvv2ff class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEvv2ff>
  : public ClassTraitsBase<Herwig::MEvv2ff> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MEvv2ff"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEvv2ff is implemented. It may also include several, space-separated,
   * libraries if the class MEvv2ff depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libHwGeneralME.so"; }
};

/// \endif

}

#include "MEvv2ff.icc"

#endif /* HERWIG_MEvv2ff_H */
