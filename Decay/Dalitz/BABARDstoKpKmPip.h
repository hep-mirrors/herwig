// -*- C++ -*-
#ifndef Herwig_BABARDstoKpKmPip_H
#define Herwig_BABARDstoKpKmPip_H
//
// This is the declaration of the BABARDstoKpKmPip class.
//

#include "WeakDalitzDecay.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The BABARDstoKpKmPip class.
 *
 * @see \ref BABARDstoKpKmPipInterfaces "The interfaces"
 * defined for BABARDstoKpKmPip.
 */
class BABARDstoKpKmPip: public WeakDalitzDecay {

public:

  /**
   * The default constructor.
   */
  BABARDstoKpKmPip();
  
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
  BABARDstoKpKmPip & operator=(const BABARDstoKpKmPip &) = delete;

private:

  /**
   *  Masses and widths of the resonances
   */
  //@{
  /**
   *  Mass of the $K^*(892)^0$
   */
  Energy mKStar_;
  
  /**
   *  Width of the $K^*(892)^0$
   */
  Energy wKStar_;
  
  /**
   *  Mass of the $\phi(102)$
   */
  Energy mPhi_;
  
  /**
   *  Width of the $\phi(102)$
   */
  Energy wPhi_;

  /**
   *  Mass of the \f$f_0(980\f$
   */
  Energy mf0_980_;

  /**
   *  Width of the \f$f_0(980\f$
   */
  Energy wf0_980_;
  
  /**
   *  Mass of the $K^*_0(1430)^0$
   */
  Energy mK0_;
  
  /**
   *  Width of the $K^*_0(1430)^0$
   */
  Energy wK0_;
  
  /**
   *  Mass of the $f_0(1710)^0$
   */
  Energy mf0_1710_;
  
  /**
   *  Width of the $f_0(1710)^0$
   */
  Energy wf0_1710_;
  
  /**
   *  Mass of the $f_0(1370)^0$
   */
  Energy mf0_1370_;
  
  /**
   *  Width of the $f_0(1370)^0$
   */
  Energy wf0_1370_;
  //@}

  /**
   *  Magnitudes and phases of the amplitudes for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$ 
   */
  //@{
  /**
   *  Amplitude of the $K^*(892)^0$
   */
  double aKStar_;
  
  /**
   *  Phase of the $K^*(892)^0$
   */
  double phiKStar_;
  
  /**
   *  Amplitude of the $\phi(102)$
   */
  double aPhi_;
  
  /**
   *  Phase of the $\phi(102)$
   */
  double phiPhi_;

  /**
   *  Amplitude of the \f$f_0(980\f$
   */
  double af0_980_;

  /**
   *  Phase of the \f$f_0(980\f$
   */
  double phif0_980_;
  
  /**
   *  Amplitude of the $K^*_0(1430)^0$
   */
  double aK0_;
  
  /**
   *  Phase of the $K^*_0(1430)^0$
   */
  double phiK0_;
  
  /**
   *  Amplitude of the $f_0(1710)^0$
   */
  double af0_1710_;
  
  /**
   *  Phase of the $f_0(1710)^0$
   */
  double phif0_1710_;
  
  /**
   *  Amplitude of the $f_0(1370)^0$
   */
  double af0_1370_;
  
  /**
   *  Phase of the $f_0(1370)^0$
   */
  double phif0_1370_;
  //@}
};

}

#endif /* Herwig_BABARDstoKpKmPip_H */
