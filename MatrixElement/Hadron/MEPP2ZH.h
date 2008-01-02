// -*- C++ -*-
#ifndef HERWIG_MEPP2ZH_H
#define HERWIG_MEPP2ZH_H
//
// This is the declaration of the MEPP2ZH class.
//

#include "Herwig++/MatrixElement/General/ME2toVH.h"
#include "MEPP2ZH.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2ZH class implements the matrix element
 * for \f$q\bar{q}\to Z^0h^0\f$.
 *
 * @see \ref MEPP2ZHInterfaces "The interfaces"
 * defined for MEPP2ZH.
 */
class MEPP2ZH: public ME2toVH {

public:

  /**
   * The default constructor.
   */
  inline MEPP2ZH();

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;
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
  static ClassDescription<MEPP2ZH> initMEPP2ZH;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2ZH & operator=(const MEPP2ZH &);

private:

  /**
   *  The allowed flavours of the incoming quarks
   */
  unsigned int _maxflavour;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2ZH. */
template <>
struct BaseClassTrait<Herwig::MEPP2ZH,1> {
  /** Typedef of the first base class of MEPP2ZH. */
  typedef Herwig::ME2toVH NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2ZH class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2ZH>
  : public ClassTraitsBase<Herwig::MEPP2ZH> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2ZH"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2ZH is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2ZH depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so"; }
};

/** @endcond */

}

#include "MEPP2ZH.icc"

#endif /* HERWIG_MEPP2ZH_H */
