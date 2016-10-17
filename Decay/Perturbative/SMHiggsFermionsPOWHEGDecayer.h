// -*- C++ -*-
#ifndef HERWIG_SMHiggsFermionsPOWHEGDecayer_H
#define HERWIG_SMHiggsFermionsPOWHEGDecayer_H
//
// This is the declaration of the SMHiggsFermionsPOWHEGDecayer class.
//

#include "SMHiggsFermionsDecayer.h"
#include "Herwig/Utilities/Maths.h"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SMHiggsFermionsPOWHEGDecayer class.
 *
 * @see \ref SMHiggsFermionsPOWHEGDecayerInterfaces "The interfaces"
 * defined for SMHiggsFermionsPOWHEGDecayer.
 */
class SMHiggsFermionsPOWHEGDecayer: public SMHiggsFermionsDecayer {

public:

  /**
   * The default constructor.
   */
  SMHiggsFermionsPOWHEGDecayer();

  /**
   *  Virtual members to be overridden by inheriting classes
   *  which implement hard corrections 
   */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return FSR;}

  /**
   *  Apply the POWHEG style correction
   */
  virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr);
  //@}

  virtual double me2(const int ichan, const Particle & part,
		     const ParticleVector & decay, MEOption meopt) const;

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SMHiggsFermionsPOWHEGDecayer> initSMHiggsFermionsPOWHEGDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHiggsFermionsPOWHEGDecayer & operator=(const SMHiggsFermionsPOWHEGDecayer &);

  /**
   *  Calcluate the Kallen function
   */
  double calculateLambda(double x, double y, double z) const;

  /**
   *  Dipole subtraction term
   */
  InvEnergy2 dipoleSubtractionTerm(double x1, double x2) const;

  /**
   *  Real emission term
   */
  InvEnergy2 calculateRealEmission(double x1, double x2) const;

  /**
   *  Virtual term
   */
  double calculateVirtualTerm() const;

  /**
   *  Non-singlet term
   */
  double calculateNonSingletTerm(double beta, double L) const;

  /**
   *  Check the sign of the momentum in the \f$z\f$-direction is correct.
   */
  bool checkZMomenta(double x1, double x2, double x3, double y, Energy pT) const;

  /**
   *  Calculate the Jacobian
   */
  InvEnergy calculateJacobian(double x1, double x2, Energy pT) const;

  /**
   *  Generate a real emission event
   */
  bool getEvent();

private:

  /**
   *  The colour factor 
   */
  double CF_;

  /**
   *  The Higgs mass
   */
  mutable Energy mHiggs_;

  /**
   *  The reduced mass
   */
  mutable double mu_;

  /**
   *  The square of the reduced mass
   */
  mutable double mu2_;

  /**
   *  The strong coupling
   */
  mutable double aS_;

  /**
   *  Stuff ofr the POWHEG correction
   */
  //@{
  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alphaS_;

  /**
   *  ParticleData object for the gluon
   */
  tcPDPtr gluon_;

  /**
   *  The cut off on pt, assuming massless quarks.
   */
  Energy pTmin_;

  //  radiative variables (pt,y)
  Energy pT_;

  /**
   *  The ParticleData objects for the fermions
   */
  vector<tcPDPtr> partons_;

  /**
   * The fermion momenta
   */
  vector<Lorentz5Momentum> quark_;

  /**
   *  The momentum of the radiated gauge boson
   */
  Lorentz5Momentum gauge_;

  /**
   *  The Higgs boson
   */
  PPtr higgs_;

  /**
   *  Higgs mass squared
   */
  Energy2 mh2_;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMHiggsFermionsPOWHEGDecayer. */
template <>
struct BaseClassTrait<Herwig::SMHiggsFermionsPOWHEGDecayer,1> {
  /** Typedef of the first base class of SMHiggsFermionsPOWHEGDecayer. */
  typedef Herwig::SMHiggsFermionsDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHiggsFermionsPOWHEGDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMHiggsFermionsPOWHEGDecayer>
  : public ClassTraitsBase<Herwig::SMHiggsFermionsPOWHEGDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMHiggsFermionsPOWHEGDecayer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SMHiggsFermionsPOWHEGDecayer is implemented. It may also include several, space-separated,
   * libraries if the class SMHiggsFermionsPOWHEGDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPerturbativeHiggsDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SMHiggsFermionsPOWHEGDecayer_H */
