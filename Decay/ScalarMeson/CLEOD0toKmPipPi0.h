// -*- C++ -*-
#ifndef Herwig_CLEOD0toKmPipPi0_H
#define Herwig_CLEOD0toKmPipPi0_H
//
// This is the declaration of the CLEOD0toKmPipPi0 class.
//

#include "WeakDalitzDecay.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The CLEOD0toKmPipPi0 class implements model of CLEO for the decay 
 * \f$D^0\to K^-\pi^+\pi^0\f$, Phys. Rev. D63 (2001) 092001.
 *
 * @see \ref CLEOD0toKmPipPi0Interfaces "The interfaces"
 * defined for CLEOD0toKmPipPi0.
 */
class CLEOD0toKmPipPi0: public WeakDalitzDecay {

public:

  /**
   * The default constructor.
   */
  CLEOD0toKmPipPi0();
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

protected:

  /**
   *  Calculate the amplitude
   */
  virtual Complex amplitude(int ichan) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CLEOD0toKmPipPi0 & operator=(const CLEOD0toKmPipPi0 &);

};

}

#endif /* Herwig_CLEOD0toKmPipPi0_H */
