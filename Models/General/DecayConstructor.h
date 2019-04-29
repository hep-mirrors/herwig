// -*- C++ -*-
//
// DecayConstructor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DecayConstructor_H
#define HERWIG_DecayConstructor_H
//
// This is the declaration of the DecayConstructor class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "NBodyDecayConstructorBase.h"
#include "Herwig/Decay/Radiation/DecayRadiationGenerator.h"
#include "DecayConstructor.fh"

namespace Herwig {
using namespace ThePEG;
  
/** 
 * The DecayConstructor class is an interfaced class that stores a 
 * vector of NBodyDecayConstructor objects and calls the appropriate 
 * function to create the decayers and decaymodes. There is also an interface
 * to add decay mode tags of the form a->b,c,...; which will not
 * be created.
 * 
 * @see \ref DecayConstructorInterfaces "The interfaces"
 * defined for DecayConstructor. 
 * @see Interfaced
 */
class DecayConstructor: public Interfaced {

public:

  /**
   * The default constructor.
   */
  DecayConstructor() : NBodyDecayConstructors_(0), 
		       _disableDMTags(0), _minBR(0.) {}

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

  /**
   * Function to create decayers
   * @param particles vector of ParticleData pointers to particles contained
   * in model
   * @param minBR minimum branching ratio for modes
   */
  void createDecayers(const vector<PDPtr> & particles, double minBR);

  /**
   * Check whether the decay mode given is one that should not be
   * created
   * @param tag The decay mode tag, a->b,c,d,...;
   */
  bool disableDecayMode(string tag) const;

  /**
   *  QED Generator
   */
  DecayRadiationGeneratorPtr QEDGenerator() {return QEDGenerator_;}

  /**
   * Vector of references to the objects that will construct the N-Body
   * decays.
   */
  const vector<NBodyDecayConstructorBasePtr> & decayConstructors() {
    return NBodyDecayConstructors_;
  }

  /**
   *  Get minimum branching ratio
   */
  double minimumBR() const { return _minBR;}

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DecayConstructor & operator=(const DecayConstructor &) = delete;

private:

  /**
   * Vector of references to the objects that will construct the N-Body
   * decays.
   */
   vector<NBodyDecayConstructorBasePtr> NBodyDecayConstructors_;

  /**
   * A list of DecayMode tags that are not to be created 
   */
  vector<string> _disableDMTags;

  /**
   *  The decay radiation generator to use for QED radiation
   */
  DecayRadiationGeneratorPtr QEDGenerator_;

  /**
   *  Minimum allowed branching ratio
   */
  double _minBR;
};

}


#endif /* HERWIG_DecayConstructor_H */
