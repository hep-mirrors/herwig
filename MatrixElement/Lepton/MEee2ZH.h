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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2ZH & operator=(const MEee2ZH &) = delete;

};

}

#endif /* HERWIG_MEee2ZH_H */
