// -*- C++ -*-
//
// MatchboxHtScale.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxHtScale_H
#define Herwig_MatchboxHtScale_H
//
// This is the declaration of the MatchboxHtScale class.
//

#include "Herwig/MatrixElement/Matchbox/Utility/MatchboxScaleChoice.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/Cuts/JetFinder.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxHtScale implements scale choices related to transverse momenta.
 *
 */
class MatchboxHtScale: public MatchboxScaleChoice {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxHtScale();

  /**
   * The destructor.
   */
  virtual ~MatchboxHtScale();
  //@}

public:

  /**
   * Return the renormalization scale. This default version returns
   * shat.
   */
  virtual Energy2 renormalizationScale() const;

  /**
   * Return the factorization scale. This default version returns
   * shat.
   */
  virtual Energy2 factorizationScale() const;

protected:

  /**
   * Initialize potential weighting factors for the given final state
   */
  virtual void initWeightFactors(const tcPDVector&, const vector<LorentzMomentum>&,
				 const Ptr<JetFinder>::ptr&) const {}

  /**
   * Return the jet pt weighting factor for the given jet
   */
  virtual double jetPtWeight(const LorentzMomentum&) const { return 1.; }

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
   * Reference to the jet finder
   */
  Ptr<JetFinder>::ptr theJetFinder;

  /**
   * Include the transverse masses of the non-jet objects
   */
  bool theIncludeMT;

  /**
   * An ovewrall scaling factor for the jet contribution
   */
  double theHTFactor;

  /**
   * An ovewrall scaling factor for the non-jet contribution
   */
  double theMTFactor;
  
  /**
   * We define jets as clustered objects with a minimum pT
   */
  Energy theScalePtCut;

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxHtScale & operator=(const MatchboxHtScale &);

};

}

#endif /* Herwig_MatchboxHtScale_H */
