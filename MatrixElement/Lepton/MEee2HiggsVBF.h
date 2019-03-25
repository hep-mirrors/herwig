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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2HiggsVBF & operator=(const MEee2HiggsVBF &) = delete;

};

}

#endif /* HERWIG_MEee2HiggsVBF_H */
