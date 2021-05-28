// -*- C++ -*-
//
// HwppSelector.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HwppSelector_H
#define HERWIG_HwppSelector_H
//
// This is the declaration of the HwppSelector class.
//

#include "HadronSelector.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup hadronization
 * The HwppSelector class selects the hadrons produced in cluster decay using
 * the Herwig variant of the cluster model.
 *
 * @see \ref HwppSelectorInterfaces "The interfaces"
 * defined for HwppSelector.
 */
class HwppSelector: public HadronSelector {

public:

  /**
   * The default constructor.
   */
  HwppSelector() : HadronSelector(1),
		   _pwtDquark( 1.0 ),_pwtUquark( 1.0 ),_pwtSquark( 1.0 ),_pwtCquark( 0.0 ),
		   _pwtBquark( 0.0 ),_pwtDIquark(1.0 ),
		   _sngWt( 1.0 ), _decWt( 1.0 ),
		   _mode(1), _enhanceSProb(0), _m0Decay(1.*GeV),
		   _scHadronWtFactor(1.), _sbHadronWtFactor(1.)
  {}

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
  pair<bool,bool> selectBaryon(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const;

  /**
   *  Strange quark weight
   */
  virtual double strangeWeight(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const;

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
  HwppSelector & operator=(const HwppSelector &) = delete;

private:

  /**
   *  The weights for the different quarks and diquarks
   */
  //@{
  /**
   * The probability of producting a down quark.
   */
  double _pwtDquark;

  /**
   * The probability of producting an up quark.
   */
  double _pwtUquark;

  /**
   * The probability of producting a strange quark.
   */
  double _pwtSquark;

  /**
   * The probability of producting a charm quark.
   */
  double _pwtCquark;

  /**
   * The probability of producting a bottom quark.
   */
  double _pwtBquark;

  /**
   * The probability of producting a diquark.
   */
  double _pwtDIquark;
  //@}

  /**
   * Singlet and Decuplet weights
   */
  //@{
  /**
   *  The singlet weight
   */
  double _sngWt;

  /**
   *  The decuplet weight
   */
  double _decWt;
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

};

}

#endif /* HERWIG_HwppSelector_H */
