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
  CLEOD0toKmPipPi0 & operator=(const CLEOD0toKmPipPi0 &) = delete;

private:

  /**
   *  Mass, Widths and related parameters
   */
  //@{

  /**
   *  Mass of the \f$K_0^*(1430)\f$
   */
  Energy mK14300_;

  /**
   *  Width of the \f$K_0^*(1430)\f$
   */
  Energy wK14300_;

  /**
   *  Mass of the \f$K^*(1680)\f$
   */
  Energy mK1680_;

  /**
   *  Width of the \f$K^*(1680)\f$
   */
  Energy wK1680_;
  /**
   *  Mass of the \f$\rho(1700)\f$
   */
  Energy mrho1700_;

  /**
   *  Width of the \f$\rho(1700)\f$
   */
  Energy wrho1700_;

  /**
   *  Mass of the \f$K^{*0}(892)\f$
   */
  Energy mK8920_;

  /**
   *  Width of the \f$K^{*0}(892)\f$
   */
  Energy wK8920_;

  /**
   *  Mass of the \f$K^{*+}(892)\f$
   */
  Energy mK892_;

  /**
   *  Width of the \f$K^{*+}(892)\f$
   */
  Energy wK892_;

  /**
   *  Mass of the \f$\rho(770)\f$
   */
  Energy mrho_;

  /**
   *  Width of the \f$\rho(770)\f$
   */
  Energy wrho_;
  //@}

  /**
   *  Magnitudes and phases of the amplitudes for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  //@{
  /**
   *  Amplitude of the non-resonant component
   */
  double aNR_;

  /**
   *  Phase of the non=resonant component
   */
  double phiNR_;

  /**
   *  Amplitude of the \f$\rho^+\f$ component
   */
  double arho_;

  /**
   *  Phase of the \f$\rho^+\f$ component
   */
  double phirho_;

  /**
   *  Amplitude of the \f$K^{*-}\f$ component
   */
  double aKstarm_;

  /**
   *  Phase of the \f$K^{*-}\f$ component
   */
  double phiKstarm_;

  /**
   *  Amplitude of the \f$\bar{K}^{*0}\f$ component
   */
  double aKstar0_;

  /**
   *  Phase of the \f$\bar{K}^{*0}\f$ component
   */
  double phiKstar0_;

  /**
   *  Amplitude for the \f$K_0(1430)^-\f$ component
   */
  Energy2 aK1430m_;

  /**
   *  Phase for the \f$K_0(1430)^-\f$ component
   */
  double phiK1430m_;

  /**
   *  Amplitude for the \f$\bar{K}_0(1430)^0\f$ component
   */
  Energy2 aK14300_;

  /**
   *  Phase for the \f$\bar{K}_0(1430)^0\f$ component
   */
  double phiK14300_;

  /**
   *  Amplitude for the \f$\rho(1700)^+\f$ component
   */
  double arho1700_;

  /**
   * Phase for the \f$\rho(1700)^+\f$ component
   */
  double phirho1700_;

  /**
   *  Amplitude of the \f$K^*(1680)^-\f$ component
   */
  double aK1680_;

  /**
   *  Phase of the \f$K^*(1680)^-\f$ component
   */
  double phiK1680_;

  /**
   *  Complex amplitude of the non-resonant component
   */
  Complex cNR_;
  //@}
  
};

}

#endif /* Herwig_CLEOD0toKmPipPi0_H */
