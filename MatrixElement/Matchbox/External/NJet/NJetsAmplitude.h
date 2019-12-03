// -*- C++ -*-
//
// NJetsAmplitude.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_NJetsAmplitude_H
#define Herwig_NJetsAmplitude_H
//
// This is the declaration of the NJetsAmplitude class.
//

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxOLPME.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief NJetsAmplitude implements an interface to NJets
 */
class NJetsAmplitude: public MatchboxOLPME {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  NJetsAmplitude();

  /**
   * The destructor.
   */
  virtual ~NJetsAmplitude();
  //@}

public:

  /**
   * Return true, if this amplitude already includes averaging over
   * incoming parton's quantum numbers.
   */
  virtual bool hasInitialAverage() const { return false; }

  /**
   * Return true, if this amplitude already includes symmetry factors
   * for identical outgoing particles.
   */
  virtual bool hasFinalStateSymmetry() const { return false; }

  /**
   * Start the one loop provider, if appropriate, giving order and
   * contract files
   */
  virtual void signOLP(const string&, const string&);

  /**
   * Start the one loop provider, if appropriate
   */
  virtual void startOLP(const string&, int& status);

  /**
   * Start the one loop provider, if appropriate. This default
   * implementation writes an BLHA 2.0 order file and starts the OLP
   */
  virtual bool startOLP(const map<pair<Process,int>,int>& procs);

  /**
   * Call OLP_EvalSubProcess and fill in the results
   */
  virtual void evalSubProcess() const;

  /**
   * Fill in results for the given colour correlator
   */
  virtual void evalColourCorrelator(pair<int,int> ij) const;

  /**
   * Return a positive helicity polarization vector for a gluon of
   * momentum p (with reference vector n) to be used when evaluating
   * spin correlations.
   */
  virtual LorentzVector<Complex> plusPolarization(const Lorentz5Momentum& p,
						  const Lorentz5Momentum& n,
						  int id = -1) const;

  /**
   * Fill in results for the given colour/spin correlator
   */
  virtual void evalSpinColourCorrelator(pair<int,int> ij) const;

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

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NJetsAmplitude & operator=(const NJetsAmplitude &) = delete;

  /**
   * Store colour correlator results
   */
  mutable vector<double> colourCorrelatorResults;

  /**
   * Store spin colour correlator results
   */
  mutable vector<double> spinColourCorrelatorResults;

  /**
   *  Location of NJETs
   */
  string NJetsPrefix_;

  /**
   *  Location of NJET librarys
   */
  string NJetsLibs_;

  /**
   * Load the NJET library
   */
  void loadNJET();

};

}

#endif /* Herwig_NJetsAmplitude_H */
