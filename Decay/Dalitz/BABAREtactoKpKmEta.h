// -*- C++ -*-
#ifndef Herwig_BABAREtactoKpKmEta_H
#define Herwig_BABAREtactoKpKmEta_H
//
// This is the declaration of the BABAREtactoKpKmEta class.
//

#include "WeakDalitzDecay.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the BABAREtactoKpKmEta class.
 *
 * @see \ref BABAREtactoKpKmEtaInterfaces "The interfaces"
 * defined for BABAREtactoKpKmEta.
 */
class BABAREtactoKpKmEta: public WeakDalitzDecay {

public:

  /**
   * The default constructor.
   */
  BABAREtactoKpKmEta();
  
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
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BABAREtactoKpKmEta & operator=(const BABAREtactoKpKmEta &) = delete;

private:

  /**
   *  Masses and widths of the resonances
   */
  //@{
  /**
   *  Mass of the $f_0(1500)^0$
   */
  Energy mf0_1500_;
  
  /**
   *  Width of the $f_0(1500)^0$
   */
  Energy wf0_1500_;
  
  /**
   *  Mass of the $f_0(1710)^0$
   */
  Energy mf0_1710_;
  
  /**
   *  Width of the $f_0(1710)^0$
   */
  Energy wf0_1710_;
  
  /**
   *  Mass of the $K^*_0(1430)^+$
   */
  Energy mK0_1430_;
  
  /**
   *  Width of the $K^*_0(1430)^+$
   */
  Energy wK0_1430_;
  
  /**
   *  Mass of the $f_0(2200)^0$
   */
  Energy mf0_2200_;
  
  /**
   *  Width of the $f_0(2200)^0$
   */
  Energy wf0_2200_;
  
  /**
   *  Mass of the $K^*_0(1950)^+$
   */
  Energy mK0_1950_;
  
  /**
   *  Width of the $K^*_0(1950)^+$
   */
  Energy wK0_1950_;
  
  /**
   *  Mass of the $f_2(1525)^0$
   */
  Energy mf2_1525_;
  
  /**
   *  Width of the $f_2(1525)^0$
   */
  Energy wf2_1525_;
  
  /**
   *  Mass of the $f_0(1370)^0$
   */
  Energy mf0_1370_;
  
  /**
   *  Width of the $f_0(1370)^0$
   */
  Energy wf0_1370_;

  /**
   *  Mass of the \f$f_0(980\f$
   */
  Energy mf0_980_;

  /**
   *  Width of the \f$f_0(980\f$
   */
  Energy wf0_980_;
  //@}

  /**
   *  Magnitudes and phases of the amplitudes
   */
  //@{
  /**
   *  Amplitude of the $f_0(1500)^0$
   */
  double af0_1500_;
  
  /**
   *  Phase of the $f_0(1500)^0$
   */
  double phif0_1500_;
  
  /**
   *  Amplitude of the $f_0(1710)^0$
   */
  double af0_1710_;
  
  /**
   *  Phase of the $f_0(1710)^0$
   */
  double phif0_1710_;
  
  /**
   *  Amplitude of the $K^*_0(1430)^+$
   */
  double aK0_1430_;
  
  /**
   *  Phase of the $K^*_0(1430)^+$
   */
  double phiK0_1430_;
  
  /**
   *  Amplitude of the $f_0(2200)^0$
   */
  double af0_2200_;
  
  /**
   *  Phase of the $f_0(2200)^0$
   */
  double phif0_2200_;
  
  /**
   *  Amplitude of the $K^*_0(1950)^+$
   */
  double aK0_1950_;
  
  /**
   *  Phase of the $K^*_0(1950)^+$
   */
  double phiK0_1950_;
  
  /**
   *  Amplitude of the $f_2(1525)^0$
   */
  double af2_1525_;
  
  /**
   *  Phase of the $f_2(1525)^0$
   */
  double phif2_1525_;
  
  /**
   *  Amplitude of the $f_0(1370)^0$
   */
  double af0_1370_;
  
  /**
   *  Phase of the $f_0(1370)^0$
   */
  double phif0_1370_;

  /**
   *  Amplitude of the \f$f_0(980\f$
   */
  double af0_980_;

  /**
   *  Phase of the \f$f_0(980\f$
   */
  double phif0_980_;

  /**
   *  Non-resonant amplitude
   */
  double aNR_;

  /**
   *  Non-resonant phase
   */
  double phiNR_;
  //@}

  /**
   *  Control over channels to check fit fractions
   */
  int channel1_, channel2_;
};

}

#endif /* Herwig_BABAREtactoKpKmEta_H */
