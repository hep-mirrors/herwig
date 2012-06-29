// -*- C++ -*-
#ifndef HERWIG_HiddenOneOneOneSplitFn_H
#define HERWIG_HiddenOneOneOneSplitFn_H
//
// This is the declaration of the HiddenOneOneOneSplitFn class.
//

#include "Herwig++/Shower/SplittingFunctions/OneOneOneSplitFn.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HiddenOneOneOneSplitFn class.
 *
 * @see \ref HiddenOneOneOneSplitFnInterfaces "The interfaces"
 * defined for HiddenOneOneOneSplitFn.
 */
class HiddenOneOneOneSplitFn: public OneOneOneSplitFn {

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
  static NoPIOClassDescription<HiddenOneOneOneSplitFn> initHiddenOneOneOneSplitFn;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HiddenOneOneOneSplitFn & operator=(const HiddenOneOneOneSplitFn &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HiddenOneOneOneSplitFn. */
template <>
struct BaseClassTrait<Herwig::HiddenOneOneOneSplitFn,1> {
  /** Typedef of the first base class of HiddenOneOneOneSplitFn. */
  typedef Herwig::OneOneOneSplitFn NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HiddenOneOneOneSplitFn class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HiddenOneOneOneSplitFn>
  : public ClassTraitsBase<Herwig::HiddenOneOneOneSplitFn> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HiddenOneOneOneSplitFn"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HiddenOneOneOneSplitFn is implemented. It may also include several, space-separated,
   * libraries if the class HiddenOneOneOneSplitFn depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HiddenOneOneOneSplitFn.so"; }
};

/** @endcond */

}

#endif /* HERWIG_HiddenOneOneOneSplitFn_H */
