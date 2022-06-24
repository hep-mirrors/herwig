// -*- C++ -*-
#ifndef Herwig_DalitzSigma_H
#define Herwig_DalitzSigma_H
//
// This is the declaration of the DalitzSigma class.
//

#include "DalitzResonance.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The DalitzSigma class implements the model of Zou and Bugg for the sigma resonance.
 */
class DalitzSigma: public DalitzResonance {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DalitzSigma()
  {}

  /**
   *  Constructor with parameters
   */
  DalitzSigma(long pid, ResonanceType::Type rtype, Energy m, Energy w,
	      unsigned int d1, unsigned int d2, unsigned int s,
	      double mag, double phi, InvEnergy rr,
	      Energy2 a, Energy b1, InvEnergy b2, Energy g4pi)
    : DalitzResonance(pid,rtype,m,w,d1,d2,s,mag,phi,rr),
      a_(a), b1_(b1),b2_(b2), g4Pi_(g4pi)
  {}

public:

  /**
   *  Return the Breit-Wigner times the form factor
   */
  virtual Complex BreitWigner(const Energy & mAB, const Energy & mA, const Energy & mB) const;

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
  DalitzSigma & operator=(const DalitzSigma &) = delete;

private :

  /**
   *  Four pion phase-space
   */
  double rho4pi(const Energy2 &s,const Energy & mpi) const {
    static const InvEnergy2 c1(3.5/GeV2);
    static const Energy2    c2(2.8*GeV2);
    return sqrt(1.-16.*sqr(mpi)/s)/(1.+exp(c1*(c2-s)));
  }
  
private:

  /**
   *  Parameters of the propagator form
   */
  //@{
  /**
   *  \f$a\f$ parameter
   */
  Energy2 a_;

  /**
   *  \f$b_1\f$ parameter
   */
  Energy b1_;

  /**
   *  \f$b_2\f$ parameter
   */
  InvEnergy b2_;

  /**
   *  Four pion coupling
   */
  Energy g4Pi_;
  //@}
  
};

}

#endif /* Herwig_DalitzSigma_H */
