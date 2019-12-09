// -*- C++ -*-
//
// MatchboxAmplitudellbarqqbargg.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxAmplitudellbarqqbargg_H
#define Herwig_MatchboxAmplitudellbarqqbargg_H
//
// This is the declaration of the MatchboxAmplitudellbarqqbargg class.
//

#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxZGammaAmplitude.h"
#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxCurrents.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxAmplitudellbarqqbargg
 */
class MatchboxAmplitudellbarqqbargg: 
    public MatchboxZGammaAmplitude, public MatchboxCurrents {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxAmplitudellbarqqbargg();

  /**
   * The destructor.
   */
  virtual ~MatchboxAmplitudellbarqqbargg();
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
  virtual unsigned int orderInGs() const { return 2; }

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const { return 2; }

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
   * Flush all cashes.
   */
  virtual void flushCaches() {
    MatchboxCurrents::reset();
    MatchboxZGammaAmplitude::flushCaches();
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
  MatchboxAmplitudellbarqqbargg & operator=(const MatchboxAmplitudellbarqqbargg &) = delete;

};

}

#endif /* Herwig_MatchboxAmplitudellbarqqbargg_H */
