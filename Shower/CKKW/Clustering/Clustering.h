// -*- C++ -*-
//
// Clustering.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Clustering_H
#define HERWIG_Clustering_H
//
// This is the declaration of the Clustering class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "Clustering.fh"

#include "ClusteringParticle.h"
#include "Clusterer.fh"
#include "PostClustering.h"
#include "ClusteringConfiguration.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * Clustering is the base class for all clusterings actually
 * being performed in a cascade history reconstruction.
 *
*@author Simon Plaetzer
 *
 * @see \ref ClusteringInterfaces "The interfaces"
 * defined for Clustering.
 */
class Clustering: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Clustering();

  /**
   * Construct a clustering giving children, data of emerging particles
   * and clusterer to perform this clustering.
   */
  inline Clustering(const vector<tClusteringParticlePtr>&,
		    const vector<ClusteringParticlePtr>&,
		    tClustererPtr,
		    tClusteringConfigurationPtr);

  /**
   * The destructor.
   */
  virtual ~Clustering();
  //@}

public:

  /**
   * Perform this clustering.
   *
   * The clustering is presented to the clusterer
   * for computing the kinematics via
   * Clusterer::doKinematics (tClusteringPtr).
   * Finally, wasClustered is called for each children and emergedFromClustering
   * for each parent.
   */
  void perform ();

  /**
   * Undo this clustering.
   *
   * Call wasUnclustered on each children particle and return
   * the particles emerging from this clustering.
   */
  vector<ClusteringParticlePtr> undo ();

public:

  /**@name Access to children and parents. */
  //@{

  /**
   * Get the children of this clustering
   */
  inline vector<tClusteringParticlePtr> children () const;

  /**
   * Get the parents of this clustering
   */
  inline vector<ClusteringParticlePtr> parents () const;

  //@}

public:

  /**@name Access to the PostClustering object */
  //@{

  /**
   * Set the PostClustering object
   */
  inline void postClustering (PostClusteringPtr);

  /**
   * Get the PostClustering object
   */
  inline tPostClusteringPtr postClustering () const;

  //@}

public:
  /**@name Methods related to clustering parameters */
  //@{

  /**
   * Set the scale of this clustering
   */
  inline void scale (const Energy2&);

  /**
   * Get the scale of this clustering
   */
  inline Energy2 scale () const;

  /**
   * Set the scale the running coupling is
   * to be evaluated for this clustering.
   */
  inline void alphaScale (const Energy2&);

  /**
   * Get the scale the running coupling is
   * to be evaluated for this clustering.
   */
  inline Energy2 alphaScale () const;

  /**
   * Set the momentumFraction of this clustering
   */
  inline void momentumFraction (double);

  /**
   * Get the momentumFraction of this clustering
   */
  inline double momentumFraction () const;

  /**
   * Set the weight of this clustering
   */
  inline void weight (double);

  /**
   * Get the weight of this clustering
   */
  inline double weight () const;

  /**
   * Veto this clustering
   */
  inline void veto ();

  /**
   * Return true, if this clustering has been vetoed
   */
  inline bool wasVetoed () const;

  /**
   * Return the clustering configuration.
   */
  inline tClusteringConfigurationPtr clusteringConfiguration () const;

  /**
   * Return the clusterer
   */
  inline tClustererPtr clusterer() const;

  //@}

public:

  /**
   * Set the event generator
   */
  inline void eventGenerator(tEGPtr);

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
   * A pointer to the clusterer which performed this clustering
   */
  tClustererPtr _clusterer;

  /**
   * The particles emerging from the clustering
   */
  vector<ClusteringParticlePtr> _parents;

  /**
   * The particles being clustered in this clustering
   */
  vector<tClusteringParticlePtr> _children;

  /**
   * The scale associated with this clustering
   */
  Energy2 _clusteringScale;

  /**
   * The scale the running coupling
   * is to be evaluated for this clustering.
   */
  Energy2 _alphaScale;

  /**
   * The momentum fraction associated with this clustering.
   */
  double _z;

  /**
   * An additional weight associated with this clustering
   */
  double _weight;

  /**
   * True, if this clustering should not be considered.
   * Set from Clusterer.
   */
  bool _veto;

  /**
   * The postclustering object to be used
   */
  PostClusteringPtr _postClustering;

  /**
   * The ClusteringConfiguration associated with
   * this clustering. This is given as a hint for
   * the clusterer and for checking consistency.
   */
  tClusteringConfigurationPtr _clusteringConfiguration;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static NoPIOClassDescription<Clustering> initClustering;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Clustering & operator=(const Clustering &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Clustering. */
template <>
struct BaseClassTrait<Herwig::Clustering,1> {
  /** Typedef of the first base class of Clustering. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Clustering class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Clustering>
  : public ClassTraitsBase<Herwig::Clustering> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::Clustering"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Clustering is implemented. It may also include several, space-separated,
   * libraries if the class Clustering depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "Clustering.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Clustering.tcc"
#endif

#endif /* HERWIG_Clustering_H */
