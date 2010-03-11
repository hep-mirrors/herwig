// -*- C++ -*-
//
// ShowerVeto.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerVeto_H
#define HERWIG_ShowerVeto_H
//
// This is the declaration of the ShowerVeto class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ShowerVeto.fh"
#include "Herwig++/Shower/SplittingFunctions/SplittingGenerator.h"
#include "Herwig++/Shower/ShowerConfig.h"
#include "ShowerTree.fh"
#include "ShowerProgenitor.fh"

namespace Herwig {

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
			     const Branching&) = 0;

  /**
   * Return true, if the selected emission off the given
   * particle and progenitor is vetoed.
   */
  virtual bool vetoSpaceLike (tcShowerProgenitorPtr, tcShowerParticlePtr,
			      const Branching&) = 0;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<ShowerVeto> initShowerVeto;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerVeto & operator=(const ShowerVeto &);

private:

  /**
   * The type of this veto.
   */
  ShowerVetoType _vetoType;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ShowerVeto. */
template <>
struct BaseClassTrait<Herwig::ShowerVeto,1> {
  /** Typedef of the first base class of ShowerVeto. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ShowerVeto class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ShowerVeto>
  : public ClassTraitsBase<Herwig::ShowerVeto> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ShowerVeto"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ShowerVeto is implemented. It may also include several, space-separated,
   * libraries if the class ShowerVeto depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ShowerVeto_H */
