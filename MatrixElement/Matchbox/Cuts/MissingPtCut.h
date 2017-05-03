// -*- C++ -*-
//
// MissingPtCut.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MissingPtCut_H
#define Herwig_MissingPtCut_H
//
// This is the declaration of the MissingPtCut class.
//

#include "ThePEG/Cuts/MultiCutBase.h"
#include "ThePEG/PDT/MatcherBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Christian Reuschle
 *
 * \brief MissingPtCut implements a cut on the total missing transverse momentum of a set of outgoing particles, i.e. for now the total transverse momentum of all outgoing neutrinos in an event.
 *
 * @see \ref MissingPtCutInterfaces "The interfaces"
 * defined for MissingPtCut.
 */
class MissingPtCut: public MultiCutBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MissingPtCut();

  /**
   * The destructor.
   */
  virtual ~MissingPtCut();
  //@}

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return the minimum allowed value of the squared invariant mass of
   * a set of outgoing partons of the given types. Typically used to
   * cut off the tails of the mass of a resonance for efficiency.
   */
  virtual Energy2 minS(const tcPDVector) const { return ZERO; }

  /**
   * Return the maximum allowed value of the squared invariant mass of
   * a set of outgoing partons of the given types. Typically used to
   * cut off the tails of the mass of a resonance for efficiency.
   */
  virtual Energy2 maxS(const tcPDVector) const { return Constants::MaxEnergy2; }

  /**
   * Return true if a set of outgoing particles with type \a ptype
   * and corresponding momenta \a p passes the cuts.
   */
  virtual bool passCuts(tcCutsPtr parent, const tcPDVector & ptype,
			const vector<LorentzMomentum> & p) const;

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;

  /**
   * Return the matcher for particles to cut on.
   */
  Ptr<MatcherBase>::tptr matcher() const { return theMatcher; }
  //@}

public:

  /**
   * Return the PDG codes of those particles that cannot be detected
   */
  const vector<int>& invisibleParticles() const { return theInvisibleParticles; }

  /**
   * Command to insert the PDG code of a particle that cannot be detected
   */
  string doInvisibleParticles(string);

  /**
   * Return the minimum missing pt.
   */
  Energy ptMissMin() const { return thePtMissMin; }

  /**
   * Return the maximum missing pt.
   */
  Energy ptMissMax() const { return thePtMissMax; }

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
   * The PDG codes of those particles that cannot be detected
   */
  vector<int> theInvisibleParticles;

  /**
   * The minimum missing pt.
   */
  Energy thePtMissMin;

  /**
   * The maximum missing pt.
   */
  Energy thePtMissMax;

  /**
   * A matcher for particles to cut on.
   */
  Ptr<MatcherBase>::ptr theMatcher;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MissingPtCut & operator=(const MissingPtCut &);

};

}

#endif /* Herwig_MissingPtCut_H */
