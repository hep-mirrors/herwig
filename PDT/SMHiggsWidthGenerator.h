// -*- C++ -*-
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
 * Here is the documentation of the SMHiggsWidthGenerator class.
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
  //@}

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

  /** @routines to calculate Higgs width. */
  //@{
  /**
   * Calculates the Higgs width with some NLL corrections a-la FORTRAN HERWIG. 
   * The following channels are taken into account: 
   * H->1quarks, H->2leptons, H->WW, H->ZZ, H->2gammas, H->2gluons
   * The prescription corresponds to one in FORTRAN HERWIG (except H->2gluons!)
   * @returns the Higgs width for the Higgs mass Mh.
   */
  Energy calcNLLRunningWidth(Energy Mh) const;

  /**
   * Calculates the Higgs width at LO.
   * @returns the Higgs width for Higgs mass Mh.
   */
  Energy calcLORunningWidth(Energy Mh) const;

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
 

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
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
   * Defines which decay modes are taken into account (see class documentation)
   */
  unsigned int _branchingopt;

  /**
   *  Storage of the partial widths
   */
  mutable vector<Energy> _partialwidths;

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
