// -*- C++ -*-
#ifndef Herwig_FlatteResonance2_H
#define Herwig_FlatteResonance2_H
//
// This is the declaration of the FlatteResonance2 class.
//

#include "DalitzResonance.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The FlatteResonance2 class implements the Flatte line-shape for resonances in Dalitz decays.
 */
class FlatteResonance2: public DalitzResonance {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FlatteResonance2() {}

  /**
   *  Constructor specifiying the parameters
   */
  FlatteResonance2(long pid, ResonanceType::Type rtype, Energy m, Energy w,
		  unsigned int d1, unsigned int d2, unsigned int s,
		  double mag, double phi, InvEnergy rr, Energy f1, Energy f2)
    : DalitzResonance(pid,rtype,m,w,d1,d2,s,mag,phi,rr),
      g1_(f1),g2_(f2)
  {}
  //@}

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
  FlatteResonance2 & operator=(const FlatteResonance2 &) = delete;

private:
   
  /**
   *  Parameters for the \f$f_0(980)\f$
   */
  //@{
  /**
   * \f$g_\pi\f$ coupling for the \f$f_0(980)\f$ width
   */
  Energy g1_;
  
  /**
   * \f$g_K\f$ coupling for the \f$f_0(980)\f$ width
   */
  Energy g2_;
  //@}

};

}

#endif /* Herwig_FlatteResonance2_H */
