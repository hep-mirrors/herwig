// -*- C++ -*-
#ifndef Herwig_CLEOD0toK0PipPim_H
#define Herwig_CLEOD0toK0PipPim_H
//
// This is the declaration of the CLEOD0toK0PipPim class.
//

#include "ScalarTo3ScalarDalitz.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The CLEOD0toK0PipPim class implements the Dalitz plot fit
 * of the CLEO collaboration for \f$D^0\to\bar{K}^0\pi^+\pi^-\f$,
 * Phys. Rev. Lett. 89 (2002) 251802
 *
 * @see \ref CLEOD0toK0PipPimInterfaces "The interfaces"
 * defined for CLEOD0toK0PipPim.
 */
class CLEOD0toK0PipPim: public ScalarTo3ScalarDalitz {

public:
  
  /**
   * The default constructor.
   */
  CLEOD0toK0PipPim();
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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
  CLEOD0toK0PipPim & operator=(const CLEOD0toK0PipPim &) = delete;

private:

  /**
   *  Mass, Widths and related parameters
   */
  //@{
  /**
   *  Mass of the \f$\rho(770)\f$
   */
  Energy mrho_;

  /**
   *  Width of the \f$\rho(770)\f$
   */
  Energy wrho_;
  
  /**
   *  Mass of the \f$\omega\f$
   */
  Energy momega_;

  /**
   *  Width of the \f$\omega\f$
   */
  Energy womega_;

  /**
   *  Mass of the \f$f_0(980)\f$
   */
  Energy mf980_;

  /**
   *  Width of the \f$f_0(980)\f$
   */
  Energy wf980_;

  /**
   * \f$g_\pi\f$ coupling for the \f$f_0(980)\f$ width
   */
  double gpi_;

  /**
   *\f$g_K\f$ coupling for the \f$f_0(980)\f$ width
   */
  double gK_;

  /**
   *  Option for handling the width of the \f$f_0(980)\f$
   */
  bool f0opt_;

  /**
   *  Mass of the \f$f_2(1270)\f$
   */
  Energy mf2_;

  /**
   *  Width of the \f$f_2(1270)\f$
   */
  Energy wf2_;

  /**
   *  Mass of the \f$f_0(1370)\f$
   */
  Energy mf1370_;

  /**
   *  Width of the \f$f_0(1370)\f$
   */
  Energy wf1370_;

  /**
   *  Mass of the \f$K^{*+}(892)\f$
   */
  Energy mK892_;

  /**
   *  Width of the \f$K^{*+}(892)\f$
   */
  Energy wK892_;
  /**
   *  Mass of the \f$K_0^*(1430)\f$
   */
  Energy mK14300_;

  /**
   *  Width of the \f$K_0^*(1430)\f$
   */
  Energy wK14300_;

  /**
   *  Mass of the \f$K_2^*(1430)\f$
   */
  Energy mK14302_;

  /**
   *  Width of the \f$K_2^*(1430)\f$
   */
  Energy wK14302_;

  /**
   *  Mass of the \f$K^*(1680)\f$
   */
  Energy mK1680_;

  /**
   *  Width of the \f$K^*(1680)\f$
   */
  Energy wK1680_;
  //@}

  /**
   *  Magnitudes and phases of the amplitudes for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$ 
   */
  //@{
  /**
   *  Amplitude for the \f$K^{*+}\f$
   */
  double aKstarp_;

  /**
   *  Phase for the \f$K^{*+}\f$
   */
  double phiKstarp_;

  /**
   *  Amplitude for the \f$\rho^0(770)\f$
   */
  double arho_;

  /**
   *  Phase for the \f$\rho^0(770)\f$
   */
  double phirho_;

  /**
   *  Amplitude for the \f$\omega\f$
   */
  double aomega_;

  /**
   *  Phase for the \f$\omega\f$
   */
  double phiomega_;

  /**
   *  Amplitude for the \f$K^{*-}\f$
   */
  double aKstarm_;

  /**
   *  Phase for the \f$K^{*-}\f$
   */
  double phiKstarm_;

  /**
   *  Amplitude for the \f$f_0(980)\f$
   */
  Energy2 af980_;

  /**
   *  Phase for the \f$f_0(980)\f$
   */
  double phif980_;

  /**
   *  Amplitude for the \f$f_2(1270)\f$
   */
  InvEnergy2 af2_;

  /**
   *  Phase for the \f$f_2(1270)\f$
   */
  double phif2_;

  /**
   *  Amplitude for the \f$f_0(1370)\f$
   */
  Energy2 af1370_;

  /**
   *  Phase for the \f$f_0(1370)\f$
   */
  double phif1370_;

  /**
   *  Amplitude for the \f$K^*_0(1430)^-\f$
   */
  Energy2 aK14300_;

  /**
   *  Phase for the \f$K^*_0(1430)^-\f$
   */
  double phiK14300_;

  /**
   *  Amplitude for the \f$K^*_2(1430)^-\f$
   */
  InvEnergy2 aK14302_;

  /**
   *  Phase for the \f$K^*_2(1430)^-\f$
   */
  double phiK14302_;

  /**
   *  Amplitude for the \f$K^*(1680)^-\f$
   */
  double aK1680_;

  /**
   *  Phase for the \f$K^*(1680)^-\f$
   */
  double phiK1680_;

  /**
   *  Amplitude of the non-resonant component
   */
  double aNR_;

  /**
   *  Phase of the non=resonant component
   */
  double phiNR_;
  
  /**
   *  Complex amplitude of the non-resonant component
   */
  Complex cNR_;
  //@}
  
};

}

#endif /* Herwig_CLEOD0toK0PipPim_H */
