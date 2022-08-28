// -*- C++ -*-
#ifndef Herwig_MIPWA_H
#define Herwig_MIPWA_H
//
// This is the declaration of the MIPWA class.
//

#include "DalitzResonance.h"
#include "Herwig/Utilities/Interpolator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MIPWA class allows the use on experimental extractions from Model Independent Partial Wave Analyses. 
 */
class MIPWA: public DalitzResonance {

public:

  /**
   * The default constructor.
   */
  MIPWA()
  {}

  /**
   *  Constructor specifiying the parameters
   */
  MIPWA(long pid, ResonanceType::Type rtype, Energy m, Energy w,
	unsigned int d1, unsigned int d2, unsigned int s,
	double mag, double phi, InvEnergy rr,
	vector<Energy> en, vector<double> mag2, vector<double> phi2)
    : DalitzResonance(pid,rtype,m,w,d1,d2,s,mag,phi,rr),
      energy_(en), mag_(mag2), phase_(phi2)
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
  MIPWA & operator=(const MIPWA &) = delete;

private:

  /**
   *  The energy points for the data table
   */
  vector<Energy> energy_;

  /**
   *  The magnitude of the amplitudes for the data table
   */
  vector<double> mag_;

  /**
   *  The phase of the amplitudes for the data table
   */
  vector<double> phase_;
  
  /**
   * Interpolators 
   */
  //@{
  /**
   *  The interpolator for the real part
   */
  mutable Interpolator<double,Energy>::Ptr iMag_;

  /**
   *  The interpolator for the imaginary part
   */
  mutable Interpolator<double,Energy>::Ptr iPhase_;
  //@}

};

}

#endif /* Herwig_MIPWA_H */
