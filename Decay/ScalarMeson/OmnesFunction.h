// -*- C++ -*-
#ifndef Herwig_OmnesFunction_H
#define Herwig_OmnesFunction_H
//
// This is the declaration of the OmnesFunction class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "OmnesFunction.fh"
#include "Herwig/Utilities/Interpolator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The OmnesFunction class provides a base class for the implementation of the Omnes function.
 *
 * @see \ref OmnesFunctionInterfaces "The interfaces"
 * defined for OmnesFunction.
 */
class OmnesFunction: public Interfaced {

public:

  /**
   *  Method to return the function value
   */
  virtual Complex D(Energy2 s) const = 0;

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  OmnesFunction & operator=(const OmnesFunction &) = delete;

};

}

#endif /* Herwig_OmnesFunction_H */
