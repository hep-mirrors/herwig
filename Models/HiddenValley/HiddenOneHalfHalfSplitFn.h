// -*- C++ -*-
#ifndef HERWIG_HiddenOneHalfHalfSplitFn_H
#define HERWIG_HiddenOneHalfHalfSplitFn_H
//
// This is the declaration of the HiddenOneHalfHalfSplitFn class.
//

#include "Herwig++/Shower/SplittingFunctions/OneHalfHalfSplitFn.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HiddenOneHalfHalfSplitFn class.
 *
 * @see \ref HiddenOneHalfHalfSplitFnInterfaces "The interfaces"
 * defined for HiddenOneHalfHalfSplitFn.
 */
class HiddenOneHalfHalfSplitFn: public OneHalfHalfSplitFn {

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
  static NoPIOClassDescription<HiddenOneHalfHalfSplitFn> initHiddenOneHalfHalfSplitFn;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HiddenOneHalfHalfSplitFn & operator=(const HiddenOneHalfHalfSplitFn &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HiddenOneHalfHalfSplitFn. */
template <>
struct BaseClassTrait<Herwig::HiddenOneHalfHalfSplitFn,1> {
  /** Typedef of the first base class of HiddenOneHalfHalfSplitFn. */
  typedef Herwig::OneHalfHalfSplitFn NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HiddenOneHalfHalfSplitFn class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HiddenOneHalfHalfSplitFn>
  : public ClassTraitsBase<Herwig::HiddenOneHalfHalfSplitFn> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HiddenOneHalfHalfSplitFn"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HiddenOneHalfHalfSplitFn is implemented. It may also include several, space-separated,
   * libraries if the class HiddenOneHalfHalfSplitFn depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HiddenOneHalfHalfSplitFn.so"; }
};

/** @endcond */

}

#endif /* HERWIG_HiddenOneHalfHalfSplitFn_H */
