// -*- C++ -*-
//
// VBFNLOAmplitude.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_VBFNLOAmplitude_H
#define Herwig_VBFNLOAmplitude_H
//
// This is the declaration of the VBFNLOAmplitude class.
//

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxOLPME.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief VBFNLOAmplitude implements an interface to VBFNLO
 */
class VBFNLOAmplitude: public MatchboxOLPME {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOAmplitude();

  /**
   * The destructor.
   */
  virtual ~VBFNLOAmplitude();
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
   * Set options
   */
  virtual bool isCS() const { return false; }
  virtual bool isExpanded() const { return true; }
  virtual bool isBDK() const { return false; }
  virtual bool isDR() const { return false; }
  virtual bool isDRbar() const { return false; }

  /**
   * Set an OLP parameter (assuming real parameters only)
   */
  void setOLPParameter(const string& name, double value) const;

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
  VBFNLOAmplitude & operator=(const VBFNLOAmplitude &);

  /**
   * Store colour correlator results
   */
  mutable vector<double> colourCorrelatorResults;

  /**
   * Store spin colour correlator results
   */
  mutable vector<double> spinColourCorrelatorResults;

};

}

#endif /* Herwig_VBFNLOAmplitude_H */
