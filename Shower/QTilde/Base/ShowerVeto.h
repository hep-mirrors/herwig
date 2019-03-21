// -*- C++ -*-
//
// ShowerVeto.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerVeto_H
#define HERWIG_ShowerVeto_H
//
// This is the declaration of the ShowerVeto class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ShowerVeto.fh"
#include "Herwig/Shower/Core/ShowerConfig.h"
#include "Herwig/Shower/Core/Base/ShowerParticle.fh"
#include "Herwig/Shower/Core/Base/ShowerProgenitor.fh"
#include "Herwig/Shower/Core/Base/ShowerTree.fh"

namespace Herwig {

struct Branching;

using namespace ThePEG;
  
/**\ingroup Shower
 * Exception class for vetoing a showering
 */
struct VetoShower { };

/**\ingroup Shower
 * ShowerVeto is a general interface for performing
 * vetoes during showering.
 *
 * @see \ref ShowerVetoInterfaces "The interfaces"
 * defined for ShowerVeto.
 */
class ShowerVeto: public Interfaced {

public:
  
  /**
   * Define types of ShowerVetoes
   */
  enum ShowerVetoType {

    /**
     * Throw away emission, if veto encountered. Set the scale to
     * the scale of vetoed emission.
     */
    Emission = 1,

    /**
     * Throw away showering
     */
    Shower,

    /**
     * Throw away event
     */
    Event
  };

public:

  /**
   * Constructor giving the behaviour of this veto
   */
  ShowerVeto (ShowerVetoType vetoType) : _vetoType(vetoType) {}

  /**
   * Return the type of this veto
   */
  ShowerVetoType vetoType () const {return _vetoType;}

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

public:

  /**
   * Return true, if the selected emission off the given
   * particle and progenitor is vetoed.
   */
  virtual bool vetoTimeLike (tcShowerProgenitorPtr, tcShowerParticlePtr,
			     const Branching&,tcShowerTreePtr) = 0;

  /**
   * Return true, if the selected emission off the given
   * particle and progenitor is vetoed.
   */
  virtual bool vetoSpaceLike (tcShowerProgenitorPtr, tcShowerParticlePtr,
			      const Branching&,tcShowerTreePtr) = 0;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerVeto & operator=(const ShowerVeto &) = delete;

private:

  /**
   * The type of this veto.
   */
  ShowerVetoType _vetoType;

};

}

#endif /* HERWIG_ShowerVeto_H */
