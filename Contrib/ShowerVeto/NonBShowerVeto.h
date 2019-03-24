// -*- C++ -*-
#ifndef Herwig_NonBShowerVeto_H
#define Herwig_NonBShowerVeto_H
//
// This is the declaration of the NonBShowerVeto class.
//

#include "Herwig/Shower/QTilde/Base/FullShowerVeto.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The NonBShowerVeto class vetos parton showers where no b (anti)quarks are produced.
 *
 * @see \ref NonBShowerVetoInterfaces "The interfaces"
 * defined for NonBShowerVeto.
 */
class NonBShowerVeto: public FullShowerVeto {

public:
  /**
   * The default constructor.
   */
  NonBShowerVeto() {}

protected:

  /**
   *  Determine whether to not to veto the shower, to be implemented in inheriting classes
   */
  virtual bool vetoShower();


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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NonBShowerVeto & operator=(const NonBShowerVeto &) = delete;

};

}

#endif /* Herwig_NonBShowerVeto_H */
