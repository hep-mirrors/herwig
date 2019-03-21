// -*- C++ -*-
#ifndef HERWIG_MEee2HiggsVBF_H
#define HERWIG_MEee2HiggsVBF_H
//
// This is the declaration of the MEee2HiggsVBF class.
//

#include "Herwig/MatrixElement/MEfftoffH.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEee2HiggsVBF class 
 *
 * @see \ref MEee2HiggsVBFInterfaces "The interfaces"
 * defined for MEee2HiggsVBF.
 */
class MEee2HiggsVBF: public MEfftoffH {

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<MEee2HiggsVBF> initMEee2HiggsVBF;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2HiggsVBF & operator=(const MEee2HiggsVBF &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEee2HiggsVBF. */
template <>
struct BaseClassTrait<Herwig::MEee2HiggsVBF,1> {
  /** Typedef of the first base class of MEee2HiggsVBF. */
  typedef Herwig::MEfftoffH NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEee2HiggsVBF class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEee2HiggsVBF>
  : public ClassTraitsBase<Herwig::MEee2HiggsVBF> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEee2HiggsVBF"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEee2HiggsVBF is implemented. It may also include several, space-separated,
   * libraries if the class MEee2HiggsVBF depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMELepton.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEee2HiggsVBF_H */
