// -*- C++ -*-
//
// PairPtCut.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_PairPtCut_H
#define Herwig_PairPtCut_H
//
// This is the declaration of the PairPtCut class.
//

#include "ThePEG/Cuts/TwoCutBase.h"
#include "ThePEG/PDT/MatcherBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Christian Reuschle
 *
 * \brief This class implements a cut on the transverse momentum of a pair of particles
 *
 * @see \ref PairPtCutInterfaces "The interfaces"
 * defined for PairPtCut.
 */
class PairPtCut: public TwoCutBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PairPtCut() : 
    theMinPt(0*GeV), theMaxPt(Constants::MaxEnergy),
    theSameFlavourOnly(false), theOppositeSignOnly(false) {}
  //@}

public:

  /** @name Overridden virtual functions defined in the base class. */
  //@{
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
   * Return the minimum allowed squared invariant mass of two outgoing
   * partons of type \a pi and \a pj.
   */
  virtual Energy2 minSij(tcPDPtr , tcPDPtr ) const { return ZERO; }

  /**
   * Return the minimum allowed value of the negative of the squared
   * invariant mass of an incoming parton of type \a pi and an
   * outgoing parton of type \a po.
   */
  virtual Energy2 minTij(tcPDPtr , tcPDPtr ) const { return ZERO; }

  /**
   * Return the minimum allowed value of \f$\Delta
   * R_{ij}=\sqrt{\Delta\eta_{ij}^2+\Delta\phi_{ij}^2}\f$ of two
   * outgoing partons of type \a pi and \a pj.
   */
  virtual double minDeltaR(tcPDPtr , tcPDPtr ) const { return ZERO; }

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
  virtual Energy minKTClus(tcPDPtr , tcPDPtr ) const { return ZERO; }

  /**
   * Return the minimum allowed value of the Durham
   * \f$k_\perp\f$-algorithms distance measure. This is defined as
   * \f$2\min(E_j^2, E_j^2)(1-\cos\theta_{ij})/\hat{s}\f$ for two
   * outgoing partons.
   */
  virtual double minDurham(tcPDPtr , tcPDPtr ) const { return ZERO; }

  //@}

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;

public:

  /**
   * Return the minimal allowed pair pT
   */
  Energy minPt() const { return theMinPt; }

  /**
   * Return the maximal allowed pair pT
   */
  Energy maxPt() const { return theMaxPt; }

  /**
   * Return whether cut acts on same-flavour fermions only
   */
  bool sameFlavourOnly() const { return theSameFlavourOnly; }

  /**
   * Return whether cut acts on opposite-sign fermions only
   */
  bool oppositeSignOnly() const { return theOppositeSignOnly; }

  /**
   * Return the matchers for a pair of particles to cut on. 
   * Only a pair of particles, matching these objects, will be affected.
   */
  Ptr<MatcherBase>::tptr firstMatcher() const { return theFirstMatcher; }
  Ptr<MatcherBase>::tptr secondMatcher() const { return theSecondMatcher; }

protected:

   /**
   * Return the family of the given PDG id number.
   */
  int family(long id) const;

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

private:

  /**
   * The minimal allowed pT cut value
   */
  Energy theMinPt;

  /**
   * The maximal allowed pT cut value
   */
  Energy theMaxPt;

  /**
   * Whether the cut is active on same-flavour fermions only 
   * (ignored for pairs not consisting of two fermions)
   */
  bool theSameFlavourOnly;

  /**
   * Whether the cut is active on opposite-sign fermions only
   * (ignored for pairs not consisting of two fermions)
   */
  bool theOppositeSignOnly;

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
  PairPtCut & operator=(const PairPtCut &);

};

}

#endif /* Herwig_PairPtCut_H */

