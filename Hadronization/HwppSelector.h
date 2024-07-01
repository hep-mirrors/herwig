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

#include "StandardModelHadronSpectrum.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup hadronization
 * The HwppSelector class selects the hadrons produced in cluster decay using
 * the Herwig variant of the cluster model.
 *
 * @see \ref HwppSelectorInterfaces "The interfaces"
 * defined for HwppSelector.
 */
class HwppSelector: public StandardModelHadronSpectrum {

public:

  /**
   * The default constructor.
   */
  HwppSelector() : StandardModelHadronSpectrum(1),
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
  std::tuple<bool,bool,bool> selectBaryon(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const;

  /**
   *  Strange quark weight
   */
  virtual double strangeWeight(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const;

  /**
   *  Returns the mass of the lightest pair of baryons.
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent
   */
  inline Energy massLightestBaryonPair(tcPDPtr ptr1, tcPDPtr ptr2) const {
    map<pair<long,long>,tcPDPair>::const_iterator lightest =
      lightestBaryons_.find(make_pair(abs(ptr1->id()),abs(ptr2->id())));
    assert(lightest!=lightestBaryons_.end());
    return lightest->second.first->mass()+lightest->second.second->mass();
  }
  
  /**
   *  Returns the lightest pair of baryons.
   * @param ptr1 is the first  constituent
   * @param ptr2 is the second constituent
   */
  tcPDPair lightestBaryonPair(tcPDPtr ptr1, tcPDPtr ptr2) const;

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
  map<pair<long,long>,tcPDPair> lightestBaryons_;
  //@}

};

}

#endif /* HERWIG_HwppSelector_H */
