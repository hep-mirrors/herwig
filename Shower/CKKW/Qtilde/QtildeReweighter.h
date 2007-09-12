// -*- C++ -*-
#ifndef HERWIG_QtildeReweighter_H
#define HERWIG_QtildeReweighter_H
//
// This is the declaration of the QtildeReweighter class.
//

#include "Herwig++/Shower/CKKW/Reweighting/DefaultReweighter.h"
#include "QtildeReweighter.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 * A reweighter for ME/PS merging with the qtilde shower.
 *
 * @author Simon Plaetzer
 *
 * @see \ref QtildeReweighterInterfaces "The interfaces"
 * defined for QtildeReweighter.
 */
class QtildeReweighter: public DefaultReweighter {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The destructor.
   */
  virtual ~QtildeReweighter();
  //@}

public:

  /**
   * Analyze the given cascade history
   * and set missing scales.
   */
  virtual void analyzeHistory (CascadeHistory);

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<QtildeReweighter> initQtildeReweighter;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QtildeReweighter & operator=(const QtildeReweighter &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QtildeReweighter. */
template <>
struct BaseClassTrait<Herwig::QtildeReweighter,1> {
  /** Typedef of the first base class of QtildeReweighter. */
  typedef Herwig::DefaultReweighter NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QtildeReweighter class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QtildeReweighter>
  : public ClassTraitsBase<Herwig::QtildeReweighter> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QtildeReweighter"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QtildeReweighter is implemented. It may also include several, space-separated,
   * libraries if the class QtildeReweighter depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "QtildeReweighter.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QtildeReweighter.tcc"
#endif

#endif /* HERWIG_QtildeReweighter_H */
