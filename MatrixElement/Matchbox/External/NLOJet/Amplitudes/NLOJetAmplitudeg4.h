// -*- C++ -*-
//
// NLOJetAmplitudeg4.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_NLOJetAmplitudeg4_H
#define Herwig_NLOJetAmplitudeg4_H
//
// This is the declaration of the NLOJetAmplitudeg4 class.
//

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetAmplitude.h"
#include "nlo++/ampg4.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer, Jan Kotanski
 *
 * \brief NLOJetAmplitudeg4
 *
 * @see \ref NLOJetAmplitudeg4Interfaces "The interfaces"
 * defined for NLOJetAmplitudeg4.
 */
class NLOJetAmplitudeg4: public Herwig::NLOJetAmplitude<0,2,0> {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  NLOJetAmplitudeg4();

  /**
   * The destructor.
   */
  virtual ~NLOJetAmplitudeg4();
  //@}

public:

  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector&) const;

  /**
   * Return the order in alpha_s
   */
  virtual unsigned int orderInAlphaS() const { return 2; }

  /**
   * Return the colour ordered subamplitude squared associated to the
   * colour structure identitfied by the given permutation of external
   * legs.
   */
  virtual double colourOrdered2(const int*, size_t) const;

  /**
   * Calculate the amplitudes for the phasespace point stored in lastXComb.
   * Call this before a derived class's action takes place.
   */
  virtual void prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr);

  /**
   * Return a ME instance appropriate for this amplitude and the given
   * subprocesses
   */
  virtual Ptr<MatchboxMEBase>::ptr makeME(const vector<PDVector>&) const;

protected:

  /**
   * Return the tree-level matrix element squared.
   * Needs to dispatch to various su3_tree(...) calls
   */
  virtual double treeLevel2(const vector<int>&) const;

  /**
   * Return the tree/oneloop interference.
   * Needs to dispatch to various su3_1loop(...) calls
   */
  virtual double treeOneLoop(const vector<int>&) const;

  /**
   * Return the colour correlated amplitudes squared.
   * Needs to dispatch to various su3_cc(...) calls
   */
  virtual double treeLevelCC(pair<int,int>,const vector<int>&) const;

  /**
   * Return the spin/colour correlated amplitude squared.
   */
  virtual pair<double,Complex> treeLevelSCC(pair<int,int>,const vector<int>&) const;

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
  NLOJetAmplitudeg4 & operator=(const NLOJetAmplitudeg4 &);

  /**
   * A pointer to the amplitude class
   */
  nlo::ampg4* theAmplitude;

};

}

#endif /* Herwig_NLOJetAmplitudeg4_H */
