// -*- C++ -*-
//
// EmitterSpectatorClustering.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_EmitterSpectatorClustering_H
#define HERWIG_EmitterSpectatorClustering_H
//
// This is the declaration of the EmitterSpectatorClustering class.
//

#include "Herwig++/Shower/CKKW/Clustering/Clustering.h"
#include "EmitterSpectatorClustering.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * EmitterSpectatorClustering is a specific implementation
 * of Clustering for algorithms working in an
 * emitter-spectator scenario.
 *
 *@author Simon Plaetzer
 *
 * @see \ref EmitterSpectatorClusteringInterfaces "The interfaces"
 * defined for EmitterSpectatorClustering.
 */
class EmitterSpectatorClustering: public Clustering {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline EmitterSpectatorClustering();

  /**
   * Construct a clustering giving children, data of emerging particles
   * and clusterer to perform this clustering. In addition,
   * specify the indices, which one of the given particles is
   * the spectator.
   */
  inline EmitterSpectatorClustering(const vector<tClusteringParticlePtr>&,
				    const unsigned int,
				    const vector<ClusteringParticlePtr>&,
				    const unsigned int,
				    tClustererPtr,
				    tClusteringConfigurationPtr);

  /**
   * The destructor.
   */
  virtual ~EmitterSpectatorClustering();
  //@}

public:

  /**
   * Return the emitter
   */
  inline tClusteringParticlePtr emitter() const;

  /**
   * Return the emission products
   */
  inline pair<tClusteringParticlePtr,tClusteringParticlePtr> emission () const;

  /**
   * Return the spectator before clustering
   */
  inline tClusteringParticlePtr spectatorBeforeClustering () const;

  /**
   * Return the spectator after clustering
   */
  inline tClusteringParticlePtr spectatorAfterClustering () const;

  /**
   * Generate a basis for sudakov clusterings.
   * The spectator is assumed to define the collinear
   * direction. The sum (difference) of the emissions
   * momenta defines the (off-shell) forward momentum.
   */
  void generateSudakovBasis ();

  /**
   * Return the basis for sudakov clusterings, (p,n).
   */
  inline pair<Lorentz5Momentum,Lorentz5Momentum> sudakovBasis () const;

public:

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

#ifdef HERWIG_DEBUG_CKKW

public:

  virtual void debugDump (ostream&);

#endif


private:

  /**
   * The index of the spectator in the childrens vector
   */
  unsigned int _childSpectator;

  /**
   * The index of the spectator in the parents vector
   */
  unsigned int _parentSpectator;

  /**
   * The last generated sudakov basis
   */
  pair<Lorentz5Momentum,Lorentz5Momentum> _sudakovBasis;  

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<EmitterSpectatorClustering> initEmitterSpectatorClustering;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EmitterSpectatorClustering & operator=(const EmitterSpectatorClustering &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of EmitterSpectatorClustering. */
template <>
struct BaseClassTrait<Herwig::EmitterSpectatorClustering,1> {
  /** Typedef of the first base class of EmitterSpectatorClustering. */
  typedef Herwig::Clustering NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EmitterSpectatorClustering class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::EmitterSpectatorClustering>
  : public ClassTraitsBase<Herwig::EmitterSpectatorClustering> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::EmitterSpectatorClustering"; }
  /**
   * The name of a file containing the dynamic library where the class
   * EmitterSpectatorClustering is implemented. It may also include several, space-separated,
   * libraries if the class EmitterSpectatorClustering depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "EmitterSpectatorClustering.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EmitterSpectatorClustering.tcc"
#endif

#endif /* HERWIG_EmitterSpectatorClustering_H */
