// -*- C++ -*-
#ifndef HERWIG_MEfv2fs_H
#define HERWIG_MEfv2fs_H
//
// This is the declaration of the MEfv2fs class.
//

#include "GeneralHardME.h"
#include "MEfv2fs.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class is designed to implement the matrix element for 
 * fermion-vector to fermion scalar. It inherits from GeneralHardME 
 * and implements the required virtual functions.
 *
 * @see @see \ref MEfv2fsInterfaces "The Interfaces"
 * defined for MEfv2fs.
 * @see GeneralHardME
 */
class MEfv2fs: public GeneralHardME {

public:

  /**
   * The default constructor.
   */
  inline MEfv2fs();

public:

  /** @name Virtual functions required by the GeneralHardME class. */
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
  static ClassDescription<MEfv2fs> initMEfv2fs;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEfv2fs & operator=(const MEfv2fs &);
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/// \if TRAITSPECIALIZATIONS

/** This template specialization informs ThePEG about the
 *  base classes of MEfv2fs. */
template <>
struct BaseClassTrait<Herwig::MEfv2fs,1> {
  /** Typedef of the first base class of MEfv2fs. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEfv2fs class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEfv2fs>
  : public ClassTraitsBase<Herwig::MEfv2fs> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MEfv2fs"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEfv2fs is implemented. It may also include several, space-separated,
   * libraries if the class MEfv2fs depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libGeneralHardME.so"; }
};

/// \endif

}

#include "MEfv2fs.icc"

#endif /* HERWIG_MEfv2fs_H */
