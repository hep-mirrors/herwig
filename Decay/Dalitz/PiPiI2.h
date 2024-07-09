// -*- C++ -*-
#ifndef Herwig_PiPiI2_H
#define Herwig_PiPiI2_H
//
// This is the declaration of the PiPiI2 class.
//

#include "DalitzResonance.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The PiPiI2 class provides an implementation of the \f$I=2\f$ s-eave for \f$\pi\pi\f$.
 */
class PiPiI2: public DalitzResonance {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PiPiI2() {};

  /**
   *  Constructor specifiying the parameters
   */
  PiPiI2(long pid, ResonanceType::Type rtype, Energy m, Energy w,
	 unsigned int d1, unsigned int d2, unsigned int s,
	 double mag, double phi, InvEnergy rr,
	 InvEnergy a, InvEnergy2 b, InvEnergy4 c, InvEnergy6 d,
	 Energy mmin, Energy mmax, double deta)
    : DalitzResonance(pid,rtype,m,w,d1,d2,s,mag,phi,rr),
      a_(a),b_(b),c_(c),d_(d),
      mmin_(mmin), mmax_(mmax), deltaEta_(deta)
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
  PiPiI2 & operator=(const PiPiI2 &) = delete;

private :

  /**
   *  Parameters for \f$\delta^2_0\f$
   */
  //@{
  /**
   *  Parameter \f$a\f$
   */
  InvEnergy a_;
  
  /**
   *  Parameter \f$b\f$
   */
  InvEnergy2 b_;
  
  /**
   *  Parameter \f$c\f$
   */
  InvEnergy4 c_;
  
  /**
   *  Parameter \f$d\f$
   */
  InvEnergy6 d_;
  //@}

  /**
   *  Parameters for \f$\eta^0\f$
   */
  //@{
  /**
   *  Minimum mass
   */
  Energy mmin_;
  
  /**
   *  Maximum mass
   */
  Energy mmax_;

  /**
   *   \f$\Delta\eta\f$ 
   */
  double deltaEta_;
  //@}
};

}

#endif /* Herwig_PiPiI2_H */
