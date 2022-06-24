// -*- C++ -*-
#ifndef Herwig_DalitzResonance_H
#define Herwig_DalitzResonance_H
//
// This is the declaration of the DalitzResonance class.
//

#include "ThePEG/Config/ThePEG.h"
#include "DalitzResonance.fh"
#include <cassert>
namespace Herwig {

using namespace ThePEG;

namespace ResonanceType {

/**
 *  Enum for the type of resonace
 */
enum Type {NonResonant=0,
	   Spin0=1,Spin1=3,Spin2=5,
	   Spin0E691=11,Spin1E691=13,Spin2E691=15,
	   BABARf0=21, Spin0Gauss=31, Flattef0=41, Spin0Complex=51, Flattea0=61, FlatteKstar0=71,
	   Sigma=81,
	   Spin1GS=23,
	   Spin0NonResonant=101, Spin1NonResonant=103, Spin2NonResonant=105,
	   Spin0MIPWA=-1, PiPiI2=-11, KMatrix=-21, LASS=-31};
}

/**
 * The DalitzResonance class provides a container class for
 * information on resonances in multi-body dalitz decays.
 */
class DalitzResonance: public Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DalitzResonance() {}

  /**
   *  Constructor specifiying the parameters
   */
  DalitzResonance(long pid, ResonanceType::Type rtype, Energy m, Energy w,
		  unsigned int d1, unsigned int d2, unsigned int s,
		  double mag, double phi, InvEnergy rr)
    : id(pid), type(rtype), mass(m),width(w),
      daughter1(d1),daughter2(d2),spectator(s),
      amp(mag*exp(Complex(0,phi))), R(rr)
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
  DalitzResonance & operator=(const DalitzResonance &) = delete;

public:

  /**
   *  PID of resonant particle
   */
  long id;

  /**
   *  Type of the resonance
   */
  ResonanceType::Type type;

  /**
   *  Mass of the resonance
   */
  Energy mass;

  /**
   *  Width of the resonance
   */
  Energy width;

  /**
   *  The children
   */
  unsigned int daughter1,daughter2;

  /**
   *   The spectactor
   */
  unsigned int spectator;

  /**
   *  The amplitude
   */
  Complex amp;

  /**
   *  Radius for the Ballt-Weisskopf formfactor
   */
  InvEnergy R;
  
};

}

#endif /* Herwig_DalitzResonance_H */
