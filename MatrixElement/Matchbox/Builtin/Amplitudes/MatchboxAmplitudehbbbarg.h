// -*- C++ -*-
//
// MatchboxAmplitudehbbbarg.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxAmplitudehbbbarg_H
#define Herwig_MatchboxAmplitudehbbbarg_H
//
// This is the declaration of the MatchboxAmplitudehbbbarg class.
//

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"
#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxCurrents.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * 
 *
 * \brief MatchboxAmplitudehbbbarg
 */
class MatchboxAmplitudehbbbarg: 
    public MatchboxAmplitude, public MatchboxCurrents {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxAmplitudehbbbarg();

  /**
   * The destructor.
   */
  virtual ~MatchboxAmplitudehbbbarg();
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
  virtual unsigned int orderInGs() const { return 1; }

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const { return 1; }

  /**
   * Return true, if this amplitude is capable of calculating one-loop
   * (QCD) corrections.
   */
  virtual bool haveOneLoop() const { return true; }

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
   * Evaluate the amplitude for the given colour tensor id and
   * helicity assignment
   */
  //virtual Complex evaluateOneLoop(size_t, const vector<int>&);

  /**
   * Return true, if one loop corrections are given in the conventions
   * of BDK.
   */
  virtual bool isBDK() const { return true; }

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

  /*
   * The interfaced up quark mass
   */
  Energy interfaceUMass; 
  /*
   * The interfaced down quark mass
   */
  Energy interfaceDMass; 
  /* 
   * The interfaced strange quark mass
   */
  Energy interfaceSMass; 
  /*
   * The interfaced charm quark mass
   */
  Energy interfaceCMass; 
  /*
   * The interfaced bottom quark mass
   */
  Energy interfaceBMass; 
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxAmplitudehbbbarg & operator=(const MatchboxAmplitudehbbbarg &);

};

}

#endif /* Herwig_MatchboxAmplitudehbbbarg_H */
