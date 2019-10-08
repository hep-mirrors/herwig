// -*- C++ -*-
#ifndef HERWIG_MEPP2ZH_H
#define HERWIG_MEPP2ZH_H
//
// This is the declaration of the MEPP2ZH class.
//

#include "Herwig/MatrixElement/MEfftoVH.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2ZH class implements the matrix element
 * for \f$q\bar{q}\to Z^0h^0\f$.
 *
 * @see \ref MEPP2ZHInterfaces "The interfaces"
 * defined for MEPP2ZH.
 */
class MEPP2ZH: public MEfftoVH {

public:

  /**
   * The default constructor.
   */
  MEPP2ZH();

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
  virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this); }
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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2ZH & operator=(const MEPP2ZH &) = delete;

};

}

#endif /* HERWIG_MEPP2ZH_H */
