// -*- C++ -*-
#ifndef HERWIG_MEee2ZH_H
#define HERWIG_MEee2ZH_H
//
// This is the declaration of the MEee2ZH class.
//

#include "Herwig/MatrixElement/MEfftoVH.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEee2ZH class implements the matrix element
 * for \f$e^+e^-\to Z^0h^0\f$.
 *
 * @see \ref MEee2ZHInterfaces "The interfaces"
 * defined for MEee2ZH.
 */
class MEee2ZH: public MEfftoVH {

public:

  /** @name Virtual functions required by the MEfftoVH class. */
  //@{
  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   *  Has not got a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return No;}

  /**
   *  Has not got an old fashioned ME correction
   */
  virtual bool hasMECorrection() {return false;}

public:

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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
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
  static NoPIOClassDescription<MEee2ZH> initMEee2ZH;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2ZH & operator=(const MEee2ZH &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEee2ZH. */
template <>
struct BaseClassTrait<Herwig::MEee2ZH,1> {
  /** Typedef of the first base class of MEee2ZH. */
  typedef Herwig::MEfftoVH NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEee2ZH class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEee2ZH>
  : public ClassTraitsBase<Herwig::MEee2ZH> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEee2ZH"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEee2ZH is implemented. It may also include several, space-separated,
   * libraries if the class MEee2ZH depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMELepton.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEee2ZH_H */
