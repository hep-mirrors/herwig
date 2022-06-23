// -*- C++ -*-
#ifndef Herwig_DalitzLASS_H
#define Herwig_DalitzLASS_H
//
// This is the declaration of the DalitzLASS class.
//

#include "DalitzResonance.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The DalitzLASS class implements the LASS parametrization of the \f$K\pi\f$ s-wave
 */
class DalitzLASS: public DalitzResonance {

public:

  /**
   * The default constructor.
   */
  DalitzLASS() {}

  /**
   *  Constructor specifiying the parameters
   */
  DalitzLASS(long pid, ResonanceType::Type rtype, Energy m, Energy w,
	     unsigned int d1, unsigned int d2, unsigned int s,
	     double mag, double phi, InvEnergy rr, unsigned int iopt,
	     double FNR, double phiNR, double Fres, double phiRes,
	     InvEnergy aScat, InvEnergy rEff)
    : DalitzResonance(pid,rtype,m,w,d1,d2,s,mag,phi,rr), opt_(iopt),
      FNR_(FNR), phiNR_(phiNR), FRes_(Fres), phiRes_(phiRes),
      aScat_(aScat), rEff_(rEff)
  {}

public:

  /**
   *  Return the Breit-Wigner times the form factor
   */
  virtual Complex BreitWigner(const Energy & mAB, const Energy & mA, const Energy & mB) const;

  /**
   *  Output the parameters
   */
  virtual void dataBaseOutput(ofstream & output);

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DalitzLASS & operator=(const DalitzLASS &) = delete;

private:

  /**
   *  Option for the form of the function
   */
  unsigned int opt_;

  /**
   *  Parameters
   */
  //@{
  /**
   *  Non-resonant magnitude
   */
  double FNR_;

  /**
   *  Non-resonant phase
   */
  double phiNR_;

  /**
   *  Resonant magnitude
   */
  double FRes_;

  /**
   *  Resonant phase
   */
  double phiRes_;

  /**
   * Scattering length
   */
  InvEnergy aScat_;

  /**
   * Effective interaction range
   */
  InvEnergy rEff_;
  //@}

};

}

#endif /* Herwig_DalitzLASS_H */
