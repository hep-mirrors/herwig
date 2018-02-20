// -*- C++ -*-
#ifndef Herwig_SudakovCutOff_H
#define Herwig_SudakovCutOff_H
//
// This is the declaration of the SudakovCutOff class.
//

#include "SudakovCutOff.fh"
#include "ThePEG/Interface/Interfaced.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The SudakovCutOff class is the base class for cut-offs in the Sudakov
 *
 * @see \ref SudakovCutOffInterfaces "The interfaces"
 * defined for SudakovCutOff.
 */
class SudakovCutOff: public Interfaced {

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SudakovCutOff & operator=(const SudakovCutOff &);

};

}

#endif /* Herwig_SudakovCutOff_H */
