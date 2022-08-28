// -*- C++ -*-
#ifndef Herwig_DalitzGS_H
#define Herwig_DalitzGS_H
//
// This is the declaration of the DalitzGS class.
//

#include "DalitzResonance.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The DalitzGS class implements the  Gounaris-Sakurai form of the propagator
 *
 */
class DalitzGS: public DalitzResonance {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DalitzGS() {}

  /**
   *  Constructor with parameters
   */
  DalitzGS(long pid, ResonanceType::Type rtype, Energy m, Energy w,
	   unsigned int d1, unsigned int d2, unsigned int s,
	   double mag, double phi, InvEnergy rr, Energy mpi);

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
  DalitzGS & operator=(const DalitzGS &) = delete;

private:
  
  /**
   *  Pion mass
   */
  Energy mpi_;

  /**
   * The function \f$\frac{\\hat{H}}{dq^2}\f$ at \f$q^2=m^2\f$ for the GS form of the
   *  Breit-Wigner
   */
  double dh_;

  /**
   * The function \f$\\hat{H}\f$ at \f$q^2=m^2\f$ for the GS form of the
   *  Breit-Wigner
   */
  Energy2 hres_;

  /**
   * The \f$H(0)\f$ parameter  for the GS form of the
   *  Breit-Wigner
   */
  Energy2 h0_;
};

}

#endif /* Herwig_DalitzGS_H */
