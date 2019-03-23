// -*- C++ -*-
#ifndef Herwig_SMTopPOWHEGDecayer_H
#define Herwig_SMTopPOWHEGDecayer_H
//
// This is the declaration of the SMTopPOWHEGDecayer class.
//

#include "SMTopDecayer.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SMTopPOWHEGDecayer class.
 *
 * @see \ref SMTopPOWHEGDecayerInterfaces "The interfaces"
 * defined for SMTopPOWHEGDecayer.
 */
class SMTopPOWHEGDecayer: public SMTopDecayer {

public:

  /**
   * The default constructor.
   */
  SMTopPOWHEGDecayer();

  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return FSR;}

  /**
   *  Apply the POWHEG style correction
   */
  virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr);

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
 /**
   *  check if event is in dead region
   */
  bool deadZoneCheck(double xw, double xg);

protected:

  /**
   *  Calculate matrix element ratio B/R
   */
  double matrixElementRatio(vector<Lorentz5Momentum> particleMomenta);

protected:

  /**
   *  Calculate momenta of t, b, W, g
   */
  bool calcMomenta(int j, Energy pT, double y, double phi, double& xg, 
		   double& xw, double& xb, double& xb_z, 
		   vector<Lorentz5Momentum>& particleMomenta);

protected:

  /**
   *  Check the calculated momenta are physical
   */
  bool psCheck(double xg, double xw);

protected:

  /**
   *  Return the momenta including the hard emission
   */
  vector<Lorentz5Momentum> hardMomenta();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMTopPOWHEGDecayer & operator=(const SMTopPOWHEGDecayer &) = delete;

private:

  /**
   *  Top quark mass
   */
  Energy mt_;

  /**
   *  Reduced \f$W^\pm\f$ mass
   */
  double w_;

  /**
   * Reduced bottom mass
   */
  double b_;

  /**
   *  Reduced \f$W^\pm\f$ mass squared
   */
  double w2_;

  /**
   * Reduced bottom mass squared
   */
  double b2_;

  /**
   *  Minimum \f$p_T\f$
   */
  Energy pTmin_;

  /**
   *  Transverse momentum of the emission
   */
  Energy pT_;

};

}

#endif /* Herwig_SMTopPOWHEGDecayer_H */
