// -*- C++ -*-
//
// MatchboxAmplitudeggttbar.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxAmplitudeggttbar_H
#define Herwig_MatchboxAmplitudeggttbar_H
//
// This is the declaration of the MatchboxAmplitudeggttbar class.
//

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"
#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxCurrents.h"
#include "MG_gg2ttx.h"
#include "HelAmps_sm.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Alix Wilcock
 *
 * \brief MatchboxAmplitudeggttbar
 */
class MatchboxAmplitudeggttbar: 
    public MatchboxAmplitude, public MatchboxCurrents {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxAmplitudeggttbar();

  /**
   * The destructor.
   */
  virtual ~MatchboxAmplitudeggttbar();
  //@}

public:

  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector&) const;

  /**
   * Return the (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGs() const {return 2; }

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const { return 0; }

  /**
   * Return true, if this amplitude will not require colour correlations.
   */
  virtual bool noCorrelations() const { return false; }  

  /**
   * Return true, if this amplitude is capable of calculating one-loop
   * (QCD) corrections.
   */
  virtual bool haveOneLoop() const { return false; }

  /**
   * Calculate the tree level amplitudes for the phasespace point
   * stored in lastXComb.
   */
  virtual void prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr);

  /**
   * Evaluate the amplitude for the given colour tensor id and
   * helicity assignment
   */
  virtual Complex evaluate(size_t, const vector<int>&, Complex&);

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  virtual Energy2 mu2() const { return lastSHat(); }

  /**
   * Flush all cashes.
   */
  virtual void flushCaches() {
    MatchboxCurrents::reset();
    MatchboxAmplitude::flushCaches();
  }

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxAmplitudeggttbar & operator=(const MatchboxAmplitudeggttbar &) = delete;

  /**
   * Containers for external particle momenta
   */
  double mom0 [4];
  double mom1 [4];
  double mom2 [4];
  double mom3 [4];

  /**
   * Fill MGParams_ with Herwig values of parameters
   */
  void setupParams(map<string, double> & MGParams);

  /**
   * Stores parameters needed in MG_process 
   */
  map<string,double>  MGParams_;

};

}

#endif /* Herwig_MatchboxAmplitudeggttbar_H */
