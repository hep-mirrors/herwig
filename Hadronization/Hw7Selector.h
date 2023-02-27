// -*- C++ -*-
//
// Hw7Selector.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Hw7Selector_H
#define HERWIG_Hw7Selector_H
//
// This is the declaration of the Hw7Selector class.
//

#include "StandardModelHadronSpectrum.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup hadronization
 * The Hw7Selector class selects the hadrons produced in cluster decay using
 * the Herwig variant of the cluster model.
 *
 * @see \ref Hw7SelectorInterfaces "The interfaces"
 * defined for Hw7Selector.
 */
class Hw7Selector: public StandardModelHadronSpectrum {

public:

  /**
   * The default constructor.
   */
  Hw7Selector() : StandardModelHadronSpectrum(1),
		   _pwtDIquarkS0( 1.0 ),_pwtDIquarkS1( 1.0 ),
		   _mode(1), _enhanceSProb(0), _m0Decay(1.*GeV),
		   _scHadronWtFactor(1.), _sbHadronWtFactor(1.)
  {}

  /**
   * Return the particle data of the diquark (anti-diquark) made by the two
   * quarks (antiquarks) par1, par2.
   * @param par1 (anti-)quark data pointer
   * @param par2 (anti-)quark data pointer
   */
  virtual PDPtr makeDiquark(tcPDPtr par1, tcPDPtr par2);

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
  
  /**
   *  Weights for baryons
   */
  virtual double baryonWeight(long id) const;

  /**
   *  Whether to select a meson or a baryon
   */
  std::tuple<bool,bool,bool> selectBaryon(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const;

  /**
   *  Strange quark weight
   */
  virtual double strangeWeight(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const;

  /**
   *   Insert a spin\f$\frac12\f$ baryon in the table
   */
  virtual void insertOneHalf(HadronInfo a, int flav1, int flav2);

  /**
   *   Insert a spin\f$\frac32\f$ baryon in the table
   */
  virtual void insertThreeHalf(HadronInfo a, int flav1, int flav2);
  
  /**
   *  Returns the mass of the lightest pair of baryons.
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent
   * @param pspin Spin in (2S+1) of the diquark
   */
  tcPDPair lightestBaryonPair(tcPDPtr ptr1, tcPDPtr ptr2, int pspin) const;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Hw7Selector & operator=(const Hw7Selector &) = delete;

private:

  /**
   *  The weights for the diquarks
   */
  //@{
  /**
   * The probability of producting a spin-0 diquark.
   */
  double _pwtDIquarkS0;

  /**
   * The probability of producting a spin-1 diquark.
   */
  double _pwtDIquarkS1;
  //@}
  
private:

  /**
   *  Which algorithm to use
   */
  unsigned int _mode;

  /**
  *  Flag that switches between no strangeness enhancement, scaling enhancement,
  *  and exponential enhancement (in numerical order)
  */
  int _enhanceSProb;

  /**
  *  Parameter that governs the strangeness enhancement scaling
  */
  Energy _m0Decay;

  /**
  *  Flag that switches between mass measures used in strangeness enhancement:
  *  cluster mass, or the lambda measure -  ( m_{clu}^2 - (m_q + m_{qbar})^2 )
  */
  int _massMeasure;

  /**
  *  Constant variable that stops the scale in strangeness enhancement from
  *  becoming too large
  */
  const double _maxScale = 20.;

  /**
  *  Heavy strange-charm hadron wight coefficient
  */
  double _scHadronWtFactor;

  /**
  *  Heavy strange-bottom hadron wight coefficient
  */
  double _sbHadronWtFactor;

  /**
   *  Caches of lightest pairs for speed
   */
  //@{
  /**
   * Masses of lightest baryon pair
   */
  map<pair<long,long>,tcPDPair> lightestBaryonsS0_,lightestBaryonsS1_;
  //@}
};

}

#endif /* HERWIG_Hw7Selector_H */
