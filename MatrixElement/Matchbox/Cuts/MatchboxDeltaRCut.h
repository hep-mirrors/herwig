// -*- C++ -*-
//
// MatchboxDeltaRCut.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxDeltaRCut_H
#define Herwig_MatchboxDeltaRCut_H
//
// This is the declaration of the MatchboxDeltaRCut class.
//

#include "ThePEG/Cuts/TwoCutBase.h"
#include "ThePEG/PDT/MatcherBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Christian Reuschle
 *
 * \brief MatchboxDeltaRCut implements cuts related to the separation in the legoplot plane
 *
 * @see \ref MatchboxDeltaRCutInterfaces "The interfaces"
 * defined for MatchboxDeltaRCut.
 */
class MatchboxDeltaRCut: public TwoCutBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxDeltaRCut();

  /**
   * The destructor.
   */
  virtual ~MatchboxDeltaRCut();
  //@}

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return the minimum allowed value of the longitudinally invariant
   * \f$k_\perp\f$-algorithms distance measure. This is defined as
   * \f$\min(p_{\perp i}, p_{\perp
   * j})\sqrt{\Delta\eta_{ij}^2+\Delta\phi_{ij}^2}\f$ for two outgoing
   * partons, or simply \f$p_{\perp i}\f$ or \f$p_{\perp j}\f$ for a
   * single outgoing parton. Returns 0 if both partons are incoming. A
   * null pointer indicates an incoming parton, hence the type of the
   * incoming parton is irrelevant.
   */
  virtual Energy minDeltaMeasureCuts(tcPDPtr , tcPDPtr ) const { return ZERO; }
  // virtual Energy minDeltaMeasureCuts(tcPDPtr pi, tcPDPtr pj) const { return ZERO; }

  /**
   * Return the minimum allowed value of the longitudinally invariant
   * \f$k_\perp\f$-algorithms distance measure. Returns ZERO.
   */
  virtual Energy minKTClus(tcPDPtr , tcPDPtr ) const { return ZERO; }
  // virtual Energy minKTClus(tcPDPtr pi, tcPDPtr pj) const { return ZERO; }

  /**
   * Return the minimum allowed squared invariant mass of two outgoing
   * partons of type \a pi and \a pj. Returns zero.
   */
  virtual Energy2 minSij(tcPDPtr , tcPDPtr ) const { return ZERO; }
  // virtual Energy2 minSij(tcPDPtr pi, tcPDPtr pj) const { return ZERO; }

  /**
   * Return the minimum allowed value of the negative of the squared
   * invariant mass of an incoming parton of type \a pi and an
   * outgoing parton of type \a po. Returns zero.
   */
  virtual Energy2 minTij(tcPDPtr , tcPDPtr ) const { return ZERO; }
  // virtual Energy2 minTij(tcPDPtr pi, tcPDPtr po) const { return ZERO; }

  /**
   * Return the minimum allowed value of \f$\Delta
   * R_{ij}=\sqrt{\Delta\eta_{ij}^2+\Delta\phi_{ij}^2}\f$ of two
   * outgoing partons of type \a pi and \a pj. Returns zero.
   */
  virtual double minDeltaR(tcPDPtr , tcPDPtr ) const { return ZERO; }
  // virtual double minDeltaR(tcPDPtr pi, tcPDPtr pj) const { return ZERO; }

  /**
   * Return the minimum allowed value of the Durham
   * \f$k_\perp\f$-algorithms distance measure. This is defined as
   * \f$2\min(E_j^2, E_j^2)(1-\cos\theta_{ij})/\hat{s}\f$ for two
   * outgoing partons. Returns zero.
   */
  virtual double minDurham(tcPDPtr , tcPDPtr ) const { return 0.0; }
  // virtual double minDurham(tcPDPtr pi, tcPDPtr pj) const { return 0.0; }

  /**
   * Return true if a pair of particles with type \a pitype and \a
   * pjtype and momenta \a pi and \a pj respectively passes the
   * cuts. \a inci and \a inj indicates if the corresponding particles
   * are incoming.
   */
  virtual bool passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
			LorentzMomentum pi, LorentzMomentum pj,
			bool inci = false, bool incj = false) const;

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;
  //@}

public:

  /**
   * Return the minimum and maximum allowed legoplot separation
   */
  double deltaRMin() const { return theDeltaRMin; }
  double deltaRMax() const { return theDeltaRMax; }

  /**
   * Return the minimum and maximum allowed rapidity separation
   */
  double deltaYMin() const { return theDeltaYMin; }
  double deltaYMax() const { return theDeltaYMax; }

  /**
   * Return the minimum and maximum allowed azimuthal separation
   */
  double deltaPhiMin() const { return theDeltaPhiMin; }
  double deltaPhiMax() const { return theDeltaPhiMax; }

  /**
   * Return the matchers for a pair of particles to cut on. Only a pair
   * of particles, matching these objects, will be affected.
   */
  Ptr<MatcherBase>::tptr firstMatcher() const { return theFirstMatcher; }
  Ptr<MatcherBase>::tptr secondMatcher() const { return theSecondMatcher; }

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The minimum and maximum allowed legoplot separation
   */
  double theDeltaRMin;
  double theDeltaRMax;

  /**
   * The minimum and maximum allowed rapidity separation
   */
  double theDeltaYMin;
  double theDeltaYMax;

  /**
   * The minimum and maximum allowed azimuthal separation
   */
  double theDeltaPhiMin;
  double theDeltaPhiMax;

  /**
   * Matchers for a pair of particles to cut on. Only a pair
   * of particles, matching these objects, will be affected.
   */
  Ptr<MatcherBase>::ptr theFirstMatcher;
  Ptr<MatcherBase>::ptr theSecondMatcher;
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxDeltaRCut & operator=(const MatchboxDeltaRCut &) = delete;

};

}

#endif /* Herwig_MatchboxDeltaRCut_H */
