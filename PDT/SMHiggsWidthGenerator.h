// -*- C++ -*-
//
// SMHiggsWidthGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMHiggsWidthGenerator_H
#define HERWIG_SMHiggsWidthGenerator_H
//
// This is the declaration of the SMHiggsWidthGenerator class.
//

#include "GenericWidthGenerator.h"

namespace Herwig {
using namespace ThePEG;

/**
 * The SMHiggsWidthGenerator class calculates the width for the Standard Model 
 * Higgs boson
 *
 * @see \ref SMHiggsWidthGeneratorInterfaces "The interfaces"
 * defined for SMHiggsWidthGenerator.
 */
class SMHiggsWidthGenerator: public GenericWidthGenerator {

public:

  /**
   * The default constructor.
   */
  SMHiggsWidthGenerator() 
    : widthopt_(2), offshell_(10.), mw_(ZERO), mz_(ZERO), gamw_(ZERO),
      gamz_(ZERO), qmass_(7,ZERO), lmass_(3,ZERO),
      sw2_(0.), ca_(0.), cf_(0.), qlast_(ZERO)
  {}

private:


  /** @name Virtual functions to be overridden from based class */
  //@{
  /**
   * Given a particle type and a mass of an instance of that particle
   * type, calculate a width.
   */
  virtual Energy width(const ParticleData &, Energy m) const;

  /**
   * Return decay map for the given particle type.
   */
  virtual DecayMap rate(const ParticleData &) const;

  /**
   * Return a decay map for a given Particle instance.
   */
  virtual DecayMap rate(const Particle &);
  //@}

  /**
   *  Return the total width and the sum of the partial widths for
   *  modes which are used
   */
  virtual pair<Energy,Energy> width(Energy, const ParticleData &) const;

  /**
   *  Calculate the partial width for a given mode
   * @param mH The Higgs masses
   * @param imode The decay mode
   */
  Energy partialWidth(int imode, Energy mH) const;

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

  /** @name Routines to calculate Higgs width. */
  //@{
  /**
   * Calculates the double Breit-Wigner Integral a-la FORTRAN HERWIG.
   * It is used in NLL corrected Higgs width for H->WW/ZZ,
   * x = (M_V/M_H)^2, y=M_V*G_V/(M_H)^2, where M_V/G_V - V-boson mass/width
   * @return the integral value.
   */
  double HwDoubleBW(double x, double y) const;

  /**
   * Calculate a loop function for the triangle vertex GGH/AAH
   * @return the loop function value as a complex number
   */
  Complex HwW2(double tau) const;
  //@}

  /**
   *  Return the branching ratios at the given scale
   * @param scale The off shell-mass of the Higgs 
   */
  DecayMap branching(Energy scale, const ParticleData &) const;

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
  static ClassDescription<SMHiggsWidthGenerator> initSMHiggsWidthGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHiggsWidthGenerator & operator=(const SMHiggsWidthGenerator &);

private:

  /**
   * Type of the Higgs width used (options: fixed, LO running, NLL corrected running)
   */
  unsigned int widthopt_;

  /**
   *  Number of times the width the Higgs is allowed to be off-shell
   */
  double offshell_;

  /**
   *  Particle properties extracted at initialization
   */
  //@{
  /**
   *  Mass of the W boson
   */
  Energy mw_;

  /**
   *  Mass of the Z boson
   */
  Energy mz_;

  /**
   *  Width of the W boson
   */
  Energy gamw_;

  /**
   *  Width of the Z boson
   */
  Energy gamz_;

  /**
   *  Masses of the quarks
   */
  vector<Energy> qmass_;

  /**
   *  Masses of the leptons
   */
  vector<Energy> lmass_;

  /**
   *  \f$\sin^2\theta_W\f$
   */
  double sw2_;

  /**
   *  The \f$C_A\f$ colour factor
   */
  double ca_;

  /**
   *  The \f$C_F\f$ colour factor
   */
  double cf_;
  //@}

  /**
   * Storage of parameters for speed
   */
  //@{
  /**
   *  The last scale
   */
  mutable Energy qlast_;

  /**
   *  \f$\Lambda_{\rm QCD}\f$
   */
  mutable Energy lambdaQCD_;

  /**
   *  The electromagnetic coupling
   */
  mutable double alphaEM_;

  /**
   *  The strong coupling
   */
  mutable double alphaS_;

  /**
   *  Correction factor for \f$H\to q\bar{q}\f$
   */
  mutable double cd_;

  /**
   *  Fermi constant
   */
  mutable Energy2 gfermiinv_;

  /**
   *  The anomalous dimension for  \f$H\to q\bar{q}\f$
   */
  mutable double gam0_;

  /**
   *  The \f$\beta\f$ function coefficient for  \f$H\to q\bar{q}\f$
   */
  mutable double beta0_;
  //@}

  /**
   *  Map between location in decay modes vector and code here
   */
  map<int,int> locMap_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMHiggsWidthGenerator. */
template <>
struct BaseClassTrait<Herwig::SMHiggsWidthGenerator,1> {
  /** Typedef of the first base class of SMHiggsWidthGenerator. */
  typedef GenericWidthGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHiggsWidthGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMHiggsWidthGenerator>
  : public ClassTraitsBase<Herwig::SMHiggsWidthGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMHiggsWidthGenerator"; }
};

/** @endcond */

}

#endif /* HERWIG_SMHiggsWidthGenerator_H */
