// -*- C++ -*-
#ifndef HERWIG_FortranModel_H
#define HERWIG_FortranModel_H
//
// This is the declaration of the FortranModel class.
//

#include "Herwig++/Shower/Base/ShowerModel.h"
#include "FortranModel.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The FortranModel class inherits from the ShowerModel class and implements the
 * checkConsistency member for the FORTRAN Herwig shower alogorithm.
 *
 * @see \ref FortranModelInterfaces "The interfaces"
 * defined for FortranModel.
 */
class FortranModel: public ShowerModel {

public:

  /**
   * The default constructor.
   */
  inline FortranModel();

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   *  The implementation of the virtual member from the base class to
   *  check that the correct objects are loaded
   */
  virtual void checkConsistency() throw(InitException);

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
  static NoPIOClassDescription<FortranModel> initFortranModel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FortranModel & operator=(const FortranModel &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FortranModel. */
template <>
struct BaseClassTrait<Herwig::FortranModel,1> {
  /** Typedef of the first base class of FortranModel. */
  typedef Herwig::ShowerModel NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FortranModel class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FortranModel>
  : public ClassTraitsBase<Herwig::FortranModel> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::FortranModel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FortranModel is implemented. It may also include several, space-separated,
   * libraries if the class FortranModel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so HwFortranShower.so"; }
};

/** @endcond */

}

#include "FortranModel.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FortranModel.tcc"
#endif

#endif /* HERWIG_FortranModel_H */
