// -*- C++ -*-
#ifndef HERWIG_RPV_H
#define HERWIG_RPV_H
//
// This is the declaration of the RPV class.
//

#include "Herwig++/Models/Susy/MSSM.h"
#include "RPV.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the RPV class.
 *
 * @see \ref RPVInterfaces "The interfaces"
 * defined for RPV.
 */
class RPV: public MSSM {

public:

  /**
   * The default constructor.
   */
  inline RPV();

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
  static ClassDescription<RPV> initRPV;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RPV & operator=(const RPV &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of RPV. */
template <>
struct BaseClassTrait<Herwig::RPV,1> {
  /** Typedef of the first base class of RPV. */
  typedef Herwig::MSSM NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the RPV class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::RPV>
  : public ClassTraitsBase<Herwig::RPV> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::RPV"; }
  /**
   * The name of a file containing the dynamic library where the class
   * RPV is implemented. It may also include several, space-separated,
   * libraries if the class RPV depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "RPV.so"; }
};

/** @endcond */

}

#include "RPV.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "RPV.tcc"
#endif

#endif /* HERWIG_RPV_H */
