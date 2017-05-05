// -*- C++ -*-
//
// NBodyDecayConstructorBase.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_NBodyDecayConstructorBase_H
#define HERWIG_NBodyDecayConstructorBase_H
//
// This is the declaration of the NBodyDecayConstructorBase class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/PDT/ParticleData.h"
#include "NBodyDecayConstructorBase.fh"
#include "PrototypeVertex.h"
#include "DecayConstructor.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * This is the base class for NBodyDecayConstructors. An N-body 
 * decay constructor should inherit from this and implement the 
 * DecayList virtual funtcion to create the decays and decayers.  
 *
 * @see \ref NBodyDecayConstructorBaseInterfaces "The interfaces"
 * defined for NBodyDecayConstructor. 
 */
class NBodyDecayConstructorBase: public Interfaced {

public:

  /**
   * The default constructor.
   */
  NBodyDecayConstructorBase() : 
    init_(true),iteration_(1), points_(1000), info_(false), 
    createModes_(true), removeOnShell_(1), excludeEffective_(true), 
    minReleaseFraction_(1e-3), maxBoson_(1), maxList_(1),
    includeTopOnShell_(false ), removeFlavourChangingVertices_(false),
    removeSmallVertices_(false), minVertexNorm_(1e-8) 
  {}

  /**
   * Function used to determine allowed decaymodes, to be implemented
   * in derived class.
   * @param particles vector of ParticleData pointers containing 
   * particles in model
   */
  virtual void DecayList(const set<PDPtr> & particles);

  /**
   * Number of outgoing lines. Required for correct ordering.
   */
  virtual unsigned int numBodies() const = 0;

  /**
   * Set the pointer to the DecayConstrcutor
   */
  void decayConstructor(tDecayConstructorPtr d) { 
    decayConstructor_ = d;
  }

  /**
   *  Remove flavour changing vertices ?
   */
  bool removeFlavourChangingVertices() const {
    return removeFlavourChangingVertices_;
  }

  /**
   *  Remove small vertices ?
   */
  bool removeSmallVertices() const {
    return removeSmallVertices_;
  }

  /**
   *  Minimum norm for vertex removal
   */
  double minVertexNorm() const {
    return minVertexNorm_;
  }

protected:
  
  /**
   *  Method to set up the decay mode, should be overidden in inheriting class
   */
  virtual void createDecayMode(vector<NBDiagram> & mode,
			       bool possibleOnShell,
			       double symfac);
  

  /**
   * Set the branching ratio of this mode. This requires 
   * calculating a new width for the decaying particle and reweighting
   * the current branching fractions.
   * @param dm The decaymode for which to set the branching ratio
   * @param pwidth The calculated width of the mode
   */
  void setBranchingRatio(tDMPtr dm, Energy pwidth);

  /**
   * Set the interfaces of the decayers depending on the flags stored.
   * @param name Fullname of the decayer in the EventGenerator
   * including the path
   */
  void setDecayerInterfaces(string name) const;

  /**
   * Whether to initialize decayers or not
   */
  bool initialize() const { return init_; }
  
  /**
   * Number of iterations if initializing (default 1)
   */
  int iteration() const { return iteration_; }

  /**
   * Number of points to do in initialization
   */
  int points() const { return points_; }

  /**
   * Whether to output information on the decayers 
   */
  bool info() const { return info_; }

  /**
   * Whether to create the DecayModes as well as the Decayer objects 
   */
  bool createDecayModes() const { return createModes_; }

  /**
   *  Maximum number of electroweak gauge bosons
   */
  unsigned int maximumGaugeBosons() const { return maxBoson_;}

  /**
   *  Maximum number of particles from the list whose decays we are calculating
   */
  unsigned int maximumList() const { return maxList_;}

  /**
   *  Minimum energy release fraction
   */ 
  double minimumReleaseFraction() const {return minReleaseFraction_;}

  /**
   * Get the pointer to the DecayConstructor object
   */
  tDecayConstructorPtr decayConstructor() const { 
    return decayConstructor_;
  }

  /**
   *  Option for on-shell particles
   */
  unsigned int removeOnShell() const { return removeOnShell_; }

  /**
   *  Check if a vertex is excluded
   */
  bool excluded(VertexBasePtr vertex) const {
    // skip an effective vertex
    if( excludeEffective_ &&
	int(vertex->orderInGs() + vertex->orderInGem()) != int(vertex->getNpoint())-2)
      return true;
    // check if explicitly forbidden
    return excludedVerticesSet_.find(vertex)!=excludedVerticesSet_.end();
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<NBodyDecayConstructorBase> 
  initNBodyDecayConstructorBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NBodyDecayConstructorBase & operator=(const NBodyDecayConstructorBase &);

private:

  /**
   * Whether to initialize decayers or not
   */
  bool init_;
  
  /**
   * Number of iterations if initializing (default 1)
   */
  int iteration_;

  /**
   * Number of points to do in initialization
   */
  int points_;

  /**
   * Whether to output information on the decayers 
   */
  bool info_;

  /**
   * Whether to create the DecayModes as well as the Decayer objects 
   */
  bool createModes_;

  /**
   *  Whether or not to remove on-shell diagrams
   */
  unsigned int removeOnShell_;

  /**
   *  Excluded Vertices
   */
  vector<VertexBasePtr> excludedVerticesVector_;

  /**
   *  Excluded Vertices
   */
  set<VertexBasePtr> excludedVerticesSet_;

  /**
   *  Excluded Particles
   */
  vector<PDPtr> excludedParticlesVector_;

  /**
   *  Excluded Particles
   */
  set<PDPtr> excludedParticlesSet_;

  /**
   *  Whether or not to exclude effective vertices
   */
  bool excludeEffective_;
  
  /**
   * A pointer to the DecayConstructor object 
   */
  tDecayConstructorPtr decayConstructor_;

  /**
   * The minimum energy release for a three-body decay as a 
   * fraction of the parent mass
   */
  double minReleaseFraction_;

  /**
   *  Maximum number of EW gauge bosons
   */
  unsigned int maxBoson_;

  /**
   *  Maximum number of particles from the decaying particle list
   */
  unsigned int maxList_;

  /**
   *  Include on-shell for \f$t\to b W\f$
   */
  bool includeTopOnShell_;

  /**
   *  Remove flavour changing vertices ?
   */
  bool removeFlavourChangingVertices_;

  /**
   *  Remove small vertices ?
   */
  bool removeSmallVertices_;

  /**
   *  Minimum norm for vertex removal
   */
  double minVertexNorm_;
};

  /** An Exception class that can be used by all inheriting classes to
   * indicate a setup problem. */
  class NBodyDecayConstructorError : public Exception {

  public:

    NBodyDecayConstructorError() : Exception() {}
    
    NBodyDecayConstructorError(const string & str,
			       Severity sev) : Exception(str,sev)
    {}
  };

}

#endif /* HERWIG_NBodyDecayConstructorBase_H */
