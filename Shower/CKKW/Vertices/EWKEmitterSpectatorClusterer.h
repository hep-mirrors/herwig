// -*- C++ -*-
//
// EWKEmitterSpectatorClusterer.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_EWKEmitterSpectatorClusterer_H
#define HERWIG_EWKEmitterSpectatorClusterer_H
//
// This is the declaration of the EWKEmitterSpectatorClusterer class.
//

#include "EWKClusterer.h"
#include "EWKEmitterSpectatorClusterer.fh"

#include "Herwig++/Shower/CKKW/Clustering/EmitterSpectatorClustering.h"
#include "Herwig++/Shower/CKKW/Clustering/EmitterSpectatorConfiguration.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 * EWK emitter spectator clusterings.
 *
 *@author Simon Plaetzer
 *
 * @see \ref EWKEmitterSpectatorClustererInterfaces "The interfaces"
 * defined for EWKEmitterSpectatorClusterer.
 */
class EWKEmitterSpectatorClusterer: public EWKClusterer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The destructor.
   */
  virtual ~EWKEmitterSpectatorClusterer();
  //@}

public:

  /**@name Methods used to determine the ClusteringGuide */
  //@{

  /**
   * Return the number of particles this clusterer clusters.
   */
  virtual inline unsigned int toClusterMultiplicity () const;

  /**
   * Return the number of particles emerging from a clustering
   * being performed by this clusterer.
   */
  virtual inline unsigned int fromClusterMultiplicity () const;

  /**
   * Return the possible clustering configurations resulting
   * from the given configuration. Return empty vector,
   * if configuration cannot be handled by this clusterer.
   */
  virtual vector<ClusteringConfigurationPtr> configurations (const vector<ClusteringParticleData>&);

  //@}

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<EWKEmitterSpectatorClusterer> initEWKEmitterSpectatorClusterer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EWKEmitterSpectatorClusterer & operator=(const EWKEmitterSpectatorClusterer &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of EWKEmitterSpectatorClusterer. */
template <>
struct BaseClassTrait<Herwig::EWKEmitterSpectatorClusterer,1> {
  /** Typedef of the first base class of EWKEmitterSpectatorClusterer. */
  typedef Herwig::EWKClusterer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EWKEmitterSpectatorClusterer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::EWKEmitterSpectatorClusterer>
  : public ClassTraitsBase<Herwig::EWKEmitterSpectatorClusterer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::EWKEmitterSpectatorClusterer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * EWKEmitterSpectatorClusterer is implemented. It may also include several, space-separated,
   * libraries if the class EWKEmitterSpectatorClusterer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "EWKEmitterSpectatorClusterer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EWKEmitterSpectatorClusterer.tcc"
#endif

#endif /* HERWIG_EWKEmitterSpectatorClusterer_H */
