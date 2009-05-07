// -*- C++ -*-
//
// NBodyDecayConstructorBase.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "DecayConstructor.fh"

namespace Herwig {

using namespace ThePEG;



/**
 *  A struct to order the particles in the same way as in the DecayMode's
 */
struct ParticleOrdering {
  /**
   *  Operator for the ordering
   * @param p1 The first ParticleData object
   * @param p2 The second ParticleData object
   */
  bool operator()(PDPtr p1, PDPtr p2) {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};

/**
 * A set of ParticleData objects ordered as for the DecayMode's
 */
typedef multiset<PDPtr,ParticleOrdering> OrderedParticles;

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
  inline NBodyDecayConstructorBase() : 
    _init(true),_iteration(1), _points(1000), _info(false), 
    _createmodes(true) {}

  /**
   * Function used to determine allowed decaymodes, to be implemented
   * in derived class.
   * @param particles vector of ParticleData pointers containing 
   * particles in model
   */
  virtual void DecayList(const vector<PDPtr> & particles) = 0;

  /**
   * Number of outgoing lines. Required for correct ordering.
   */
  virtual unsigned int numBodies() const = 0;

  /**
   * Set the pointer to the DecayConstrcutor
   */
  inline void decayConstructor(tDecayConstructorPtr d) { 
    _decayConstructor = d;
  }

protected:
  
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
  inline bool initialize() const { return _init; }
  
  /**
   * Number of iterations if initializing (default 1)
   */
  inline int iteration() const { return _iteration; }

  /**
   * Number of points to do in initialization
   */
  inline int points() const { return _points; }

  /**
   * Whether to output information on the decayers 
   */
  inline bool info() const { return _info; }

  /**
   * Whether to create the DecayModes as well as the Decayer objects 
   */
  inline bool createDecayModes() const { return _createmodes; }

  /**
   * Get the pointer to the DecayConstructor object
   */
  inline tDecayConstructorPtr decayConstructor() const { 
    return _decayConstructor;
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
  bool _init;
  
  /**
   * Number of iterations if initializing (default 1)
   */
  int _iteration;

  /**
   * Number of points to do in initialization
   */
  int _points;

  /**
   * Whether to output information on the decayers 
   */
  bool _info;

  /**
   * Whether to create the DecayModes as well as the Decayer objects 
   */
  bool _createmodes;
  
  /**
   * A pointer to the DecayConstructor object 
   */
  tDecayConstructorPtr _decayConstructor;
};

  /** An Exception class that can be used by all inheriting classes to
   * indicate a setup problem. */
  class NBodyDecayConstructorError : public Exception {};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NBodyDecayConstructorBase. */
template <>
struct BaseClassTrait<Herwig::NBodyDecayConstructorBase,1> {
  /** Typedef of the first base class of NBodyDecayConstructorBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NBodyDecayConstructorBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NBodyDecayConstructorBase>
  : public ClassTraitsBase<Herwig::NBodyDecayConstructorBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::NBodyDecayConstructorBase"; }
};

/** @endcond */

}

#endif /* HERWIG_NBodyDecayConstructorBase_H */
