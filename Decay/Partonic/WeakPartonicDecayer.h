// -*- C++ -*-
//
// WeakPartonicDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_WeakPartonicDecayer_H
#define HERWIG_WeakPartonicDecayer_H
//
// This is the declaration of the WeakPartonicDecayer class.
//

#include "PartonicDecayerBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The WeakPartonicDecayer class is designed to replace the HeavyDecayer
 * class which implements the partonic decays of hadrons containing a 
 * heavy quark in the same way as in FORTRAN Herwig. 
 *
 * There are a number of major changes
 *
 * - The helicity formalism is used for the decays so that the \f$\tau\f$ lepton
 *   gets the correct correlations.
 *
 * - The particles produced directly by the hadronisation, i.e. the primary hadrons
 *   produced in cluster decay are checked to ensure that none of the exclusive
 *   modes are reproduced.
 *
 * - Two body modes are allowed to try and force baryon production etc. In this case
 *   the colours of the partons are connected.
 * 
 * - Three body modes of the form \f$q g \bar{q}\f$ are supported for penguin mediated
 *   weak decays.
 *
 *  Two types of matrix element are supported for this decay
 *
 *  - MECode=0   flat-phase space.
 *  - MECode=100 V-A matrix element for the heavy quark decay in the spectator model.
 *
 *  In addition for the two-body decays and the three-bopdy spectator decays 
 *  using the weka V-A matrix element the option of adding an extra gluon to increase
 *  the multiplicity of hadrons is included
 * 
 * @see HeavyDecayer
 */
class WeakPartonicDecayer: public PartonicDecayerBase {

public:

  /**
   * The default constructor.
   */
  WeakPartonicDecayer();

  /**
   * Check if this decayer can perfom the decay for a particular mode
   * @param parent The decaying particle
   * @param children The decay products
   * @return true If this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;

  
  /**
   *  Perform the decay of the particle to the specified decay products
   * @param parent The decaying particle
   * @param children The decay products
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const tPDVector & children) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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

public:

  /**
   * Weighting of phase space for V-A matrix elements
   */
  static double VAWt(Energy2, Energy2, Energy2, InvEnergy4);

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

  /**
   *  Compute \f$\rho\f$ matrix for tau decay in three body case
   * @param dec ParticleData object of decaying quark
   * @param pdec The momentum of the decaying particle
   * @param partons The partons produced
   */
  void threeBodyMatrixElement(tcPDPtr dec,Lorentz5Momentum & pdec,
			      ParticleVector& partons) const;

  /**
   *  Four body matrix element for weak decay including an extra gluon
   * @param p0 Momentum of decaying quark
   * @param p1 Momentum of connected decay product
   * @param p2 Momentum of first parton from W decay
   * @param p3 Momentum of second parton from W decay
   * @param pg Momentum of gluon from W decay
   * @param Wcol Whether or not W products are coloured
   * @param initial Whether the radiation is from the decaying quark/first decay product
   * or W decay products
   */
  double fourBodyMatrixElement(Lorentz5Momentum & p0,Lorentz5Momentum & p1,
			       Lorentz5Momentum & p2,Lorentz5Momentum & p3,
			       Lorentz5Momentum & pg,bool Wcol, bool & initial) const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<WeakPartonicDecayer> initWeakPartonicDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  WeakPartonicDecayer & operator=(const WeakPartonicDecayer &);

private:

  /**
   *  The code for the matrix element being used.
   */
  int MECode;

  /**
   *  Probablilty of radiation giving an extra quark-antiquark pair
   */
  double _radprob;

  /**
   *  Maximum number of tries to generate the kinematics
   */
  unsigned int _maxtry;

  /**
   *  Maximum weight for three-body decays
   */
  double _threemax;

  /**
   *  Maximum weight for four-body decays
   */
  double _fourmax;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of WeakPartonicDecayer. */
template <>
struct BaseClassTrait<Herwig::WeakPartonicDecayer,1> {
  /** Typedef of the first base class of WeakPartonicDecayer. */
  typedef Herwig::PartonicDecayerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the WeakPartonicDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::WeakPartonicDecayer>
  : public ClassTraitsBase<Herwig::WeakPartonicDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::WeakPartonicDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the WeakPartonicDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwPartonicDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_WeakPartonicDecayer_H */
