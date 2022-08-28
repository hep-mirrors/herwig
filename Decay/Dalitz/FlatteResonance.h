// -*- C++ -*-
#ifndef Herwig_FlatteResonance_H
#define Herwig_FlatteResonance_H
//
// This is the declaration of the FlatteResonance class.
//

#include "DalitzResonance.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The FlatteResonance class implements the Flatte line-shape for resonances in Dalitz decays.
 */
class FlatteResonance: public DalitzResonance {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FlatteResonance() {}

  /**
   *  Constructor specifiying the parameters
   */
  FlatteResonance(long pid, ResonanceType::Type rtype, Energy m, Energy w,
		  unsigned int d1, unsigned int d2, unsigned int s,
		  double mag, double phi, InvEnergy rr, vector<Energy> f)
    : DalitzResonance(pid,rtype,m,w,d1,d2,s,mag,phi,rr),
      g_(f)
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
  FlatteResonance & operator=(const FlatteResonance &) = delete;

private:
   
  /**
   * Coupling for the width
   */
  vector<Energy> g_;

};

}

#endif /* Herwig_FlatteResonance_H */
