// -*- C++ -*-
#ifndef Herwig_ScalarAmplitude_H
#define Herwig_ScalarAmplitude_H
//
// This is the declaration of the ScalarAmplitude class.
//

#include "ScalarAmplitude.h"
#include "ThePEG/Interface/Interfaced.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Base class for scalar amplitudes
 *
 * @see \ref ScalarAmplitudeInterfaces "The interfaces"
 * defined for ScalarAmplitude.
 */
class ScalarAmplitude: public Interfaced {

public:

  /**
   * The default constructor.
   */
  ScalarAmplitude();

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
  ScalarAmplitude & operator=(const ScalarAmplitude &) = delete;

};

}

#endif /* Herwig_ScalarAmplitude_H */
