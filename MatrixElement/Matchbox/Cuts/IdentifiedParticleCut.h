// -*- C++ -*-
//
// IdentifiedParticleCut.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_IdentifiedParticleCut_H
#define Herwig_IdentifiedParticleCut_H
//
// This is the declaration of the IdentifiedParticleCut class.
//

#include "ThePEG/Cuts/OneCutBase.h"
#include "ThePEG/PDT/MatcherBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief IdentifiedParticleCut implements cuts on single momenta.
 *
 * @see \ref IdentifiedParticleCutInterfaces "The interfaces"
 * defined for IdentifiedParticleCut.
 */
class IdentifiedParticleCut: public OneCutBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  IdentifiedParticleCut();

  /**
   * The destructor.
   */
  virtual ~IdentifiedParticleCut();
  //@}

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return the minimum allowed value of the transverse momentum of an
   * outgoing parton.
   */
  virtual Energy minKT(tcPDPtr) const { return ZERO; }

  /**
   * Return the minimum allowed pseudo-rapidity of an outgoing parton
   * of the given type. The pseudo-rapidity is measured in the lab
   * system.
   */
  virtual double minEta(tcPDPtr) const { return -Constants::MaxRapidity; }

  /**
   * Return the maximum allowed pseudo-rapidity of an outgoing parton
   * of the given type. The pseudo-rapidity is measured in the lab
   * system.
   */
  virtual double maxEta(tcPDPtr) const { return Constants::MaxRapidity; }

  /**
   * Return true if a particle with type \a ptype and momentum \a p
   * passes the cuts. The \a parent contains information about the
   * kinematics of the hard sub-process.
   */
  virtual bool passCuts(tcCutsPtr parent,
			tcPDPtr ptype, LorentzMomentum p) const;

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;
  //@}

public:

  /**
   * Return the minimum pt.
   */
  Energy ptMin() const { return thePtMin; }

  /**
   * Return the maximum pt.
   */
  Energy ptMax() const { return thePtMax; }

  /**
   * Return the rapidity ranges.
   */
  const vector<pair<double,double> >& yRanges() const { return theYRanges; }

  /**
   * Return the matcher for particles to cut on.
   */
  Ptr<MatcherBase>::tptr matcher() const { return theMatcher; }

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
   * Command to insert a rapidity range
   */
  string doYRange(string);

  /**
   * The minimum pt.
   */
  Energy thePtMin;

  /**
   * The maximum pt.
   */
  Energy thePtMax;

  /**
   * The rapidity ranges.
   */
  vector<pair<double,double> > theYRanges;

  /**
   * A matcher for particles to cut on.
   */
  Ptr<MatcherBase>::ptr theMatcher;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  IdentifiedParticleCut & operator=(const IdentifiedParticleCut &);

};

}

#endif /* Herwig_IdentifiedParticleCut_H */
