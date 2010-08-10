// -*- C++ -*-
#ifndef HERWIG_HiddenHalfHalfOneSplitFn_H
#define HERWIG_HiddenHalfHalfOneSplitFn_H
//
// This is the declaration of the HiddenHalfHalfOneSplitFn class.
//

#include "Herwig++/Shower/SplittingFunctions/HalfHalfOneSplitFn.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HiddenHalfHalfOneSplitFn class.
 *
 * @see \ref HiddenHalfHalfOneSplitFnInterfaces "The interfaces"
 * defined for HiddenHalfHalfOneSplitFn.
 */
class HiddenHalfHalfOneSplitFn: public HalfHalfOneSplitFn {

public:

  /**
   *  Method to check the colours are correct
   */
  virtual bool checkColours(const IdList & ids) const;

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
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<HiddenHalfHalfOneSplitFn> initHiddenHalfHalfOneSplitFn;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HiddenHalfHalfOneSplitFn & operator=(const HiddenHalfHalfOneSplitFn &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HiddenHalfHalfOneSplitFn. */
template <>
struct BaseClassTrait<Herwig::HiddenHalfHalfOneSplitFn,1> {
  /** Typedef of the first base class of HiddenHalfHalfOneSplitFn. */
  typedef Herwig::HalfHalfOneSplitFn NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HiddenHalfHalfOneSplitFn class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HiddenHalfHalfOneSplitFn>
  : public ClassTraitsBase<Herwig::HiddenHalfHalfOneSplitFn> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HiddenHalfHalfOneSplitFn"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HiddenHalfHalfOneSplitFn is implemented. It may also include several, space-separated,
   * libraries if the class HiddenHalfHalfOneSplitFn depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HiddenHalfHalfOneSplitFn.so"; }
};

/** @endcond */

}

#endif /* HERWIG_HiddenHalfHalfOneSplitFn_H */
