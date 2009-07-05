// -*- C++ -*-
//
// SMHiggsWidthGenerator.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMHiggsWidthGenerator_H
#define HERWIG_SMHiggsWidthGenerator_H
//
// This is the declaration of the SMHiggsWidthGenerator class.
//

#include "ThePEG/PDT/WidthGenerator.h"
#include "SMHiggsWidthGenerator.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * Typedef to define a DecayMoap
 */
typedef Selector<tDMPtr> DecayMap;

/**
 * The SMHiggsWidthGenerator class calculates the width for the Standard Model Higgs
 * boson.
 *
 * @see \ref SMHiggsWidthGeneratorInterfaces "The interfaces"
 * defined for SMHiggsWidthGenerator.
 */
class SMHiggsWidthGenerator: public WidthGenerator {

public:

  /**
   * The default constructor.
   */
  inline SMHiggsWidthGenerator();

  /** @name Virtual functions to be overridden from based class */
  //@{
  /**
   * Return true if this object can be used for the given particle
   * type with the given decay map.
   */
  virtual bool accept(const ParticleData &) const;

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
  pair<Energy,Energy> width(Energy, const ParticleData &) const;

  /**
   *  Calculate the partial width for a given mode
   * @param mH The Higgs masses
   * @param imode The decay mode
   */
  Energy partialWidth(Energy mH,unsigned int imode) const;

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
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
  unsigned int _widthopt;

  /**
   *  Number of times the width the Higgs is allowed to be off-shell
   */
  double _offshell;

  /**
   *  Particle properties extracged at initialization
   */
  //@{
  /**
   *  Mass of the W boson
   */
  Energy _mw;

  /**
   *  Mass of the Z boson
   */
  Energy _mz;

  /**
   *  Width of the W boson
   */
  Energy _gamw;

  /**
   *  Width of the Z boson
   */
  Energy _gamz;

  /**
   *  Masses of the quarks
   */
  vector<Energy> _qmass;

  /**
   *  Masses of the leptons
   */
  vector<Energy> _lmass;

  /**
   *  \f$\sin^2\theta_W\f$
   */
  double _sw2;

  /**
   *  The \f$C_A\f$ colour factor
   */
  double _ca;

  /**
   *  The \f$C_F\f$ colour factor
   */
  double _cf;
  //@}

  /**
   * Storage of parameters for speed
   */
  //@{
  /**
   *  The last scale
   */
  mutable Energy _qlast;

  /**
   *  \f$\Lambda_{\rm QCD}\f$
   */
  mutable Energy _lambdaQCD;

  /**
   *  The electromagnetic coupling
   */
  mutable double _alphaEM;

  /**
   *  The strong coupling
   */
  mutable double _alphaS;

  /**
   *  Correction factor for \f$H\to q\bar{q}\f$
   */
  mutable double _cd;

  /**
   *  Fermi constant
   */
  mutable Energy2 _gfermiinv;

  /**
   *  The anomalous dimension for  \f$H\to q\bar{q}\f$
   */
  mutable double _gam0;

  /**
   *  The \f$\beta\f$ function coefficient for  \f$H\to q\bar{q}\f$
   */
  mutable double _beta0;
  //@}

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
  typedef WidthGenerator NthBase;
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

#include "SMHiggsWidthGenerator.icc"

#endif /* HERWIG_SMHiggsWidthGenerator_H */
