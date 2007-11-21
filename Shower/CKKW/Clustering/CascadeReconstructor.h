// -*- C++ -*-
//
// CascadeReconstructor.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_CascadeReconstructor_H
#define HERWIG_CascadeReconstructor_H
//
// This is the declaration of the CascadeReconstructor class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/EventRecord/Particle.h"
#include "CascadeReconstructor.fh"

#include "CKKWHardProcess.h"
#include "Clusterer.h"
#include "ClusteringConfiguration.h"
#include "ClusteringGuide.h"
#include "Clustering.h"
#include "ClusteringParticle.h"
#include "ClusteringSelector.h"
#include "Partitioner.h"
#include "PostClustering.h"
#include "Herwig++/Shower/CKKW/Reweighting/JetMeasure.h"

namespace Herwig {

using namespace ThePEG;

  /**\ingroup CKKW
   *
   * Struct to organize result from
   * cascade history reconstruction.
   *
   *@author Simon Plaetzer
   *
   *@see CascadeReconstructor
   */
  struct CascadeHistory {

    /**
     * The sequence of clusterings done.
     * front() is the first one done, back() the last one.
     */
    list<ClusteringPtr> clusterings;

    /**
     * All reconstructed particles.
     */
    vector<ClusteringParticlePtr> reconstructed;

    /**
     * The hard process configuration reached.
     */
    vector<tClusteringParticlePtr> hardProcess;

    /**
     * The hard process object
     */
    tCKKWHardProcessPtr hard;

#ifdef HERWIG_DEBUG_CKKW_GRAPHVIZ

    /**
     * Output the cascade history to a file
     * to be processed by dot.
     */
    void toDot (ostream&, long evtNum, double ckkwweight);

    /**
     * Map colour line indices to X11 colour
     * scheme names.
     */
    string colour (unsigned int);

#endif

  };

/**\ingroup CKKW
 *
 * CascadeReconstructor is the class which performs reconstruction
 * of parton shower histories for ME/PS merging approaches.
 *
 *@author Simon Plaetzer
 *
 * @see \ref CascadeReconstructorInterfaces "The interfaces"
 * defined for CascadeReconstructor.
 */
class CascadeReconstructor: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline CascadeReconstructor();

  /**
   * The destructor.
   */
  virtual ~CascadeReconstructor();
  //@}

public:

  /**
   * Return the clustering particle corresponding to
   * the given ThePEG::Particle
   */
  inline tClusteringParticlePtr clusteringParticle (tPPtr) const;

  /**
   * Reconstruct cascade history starting from given
   * configuration. If cascade history exists, return true.
   */
  bool reconstruct (PPair, pair<double,double>, ParticleVector);

  /**
   * If succesfull reconstruction, return the CascadeHistory
   */
  CascadeHistory history () const;

  /**
   * Access the partitioner
   */
  inline Partitioner& partitioner ();

  /**@name Methods related to ordering of clustering scales */
  //@{

  /**
   * Return true, if increasing clustering scales are to be forced.
   */
  inline bool forceIncreasing () const;

  /**
   * Return true, if a non-ordered history should be accepted
   * in case an ordered one is non-existing.
   */
  inline bool mayUseUnordered () const;

  //@}

  /**
   * Set the jet resolution used
   */
  inline void resolution (JetMeasurePtr);

  /**
   * Do some consistency checking
   */
  void setup ();

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

  void dumpThePEGParticle (ostream&,PPtr);

#endif

protected:

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}



private:

  /**
   * Convert to ClusteringParticle's and keep a
   * conversion map.
   */
  vector<ClusteringParticlePtr> convert (PPair, pair<double,double>, ParticleVector);

  /**@name Methods used for cascade reconstruction */
  //@{
  /**
   * Recursively reconstruct.
   * Return true, if hard process reached (_currentParticles is
   * the hard process configuration), false if no clusterings
   * could be performed.
   */
  bool reconstruct ();

  /**
   * Perform the next clustering
   */
  void perform ();

  /**
   * Undo the last clustering
   */
  void undo ();

  //@}

  /**
   * The partitioner used by this CascadeReconstructor
   */
  Partitioner _partitioner;

  /**
   * A map of already looked up ClusteringGuides
   */
  map<vector<ClusteringParticleData>, ClusteringGuidePtr> _guides;

  /**
   * The clusterers to be used
   */
  vector<ClustererPtr> _clusterers;

  /**
   * The hard processes to be looked at
   */
  vector<CKKWHardProcessPtr> _hardProcesses;

  /**
   * The ClusteringSelector to be used
   */
  ClusteringSelectorPtr _selector;

  /**
   * True, if increasing clustering scales
   * should be forced.
   */
  bool _forceIncreasing;

  /**
   * Use the first available non-ordered history
   * if an ordered one is not available (only relevant
   * if _forceIncreasing = true
   */
  bool _mayUseUnordered;

  /**
   * The jet resolution used
   */
  JetMeasurePtr _resolution;

  /**@name Members associated with the current reconstruction. */
  //@{

  /**
   * All particles in the reconstruction.
   */
  vector<ClusteringParticlePtr> _particles;

  /**
   * The list of clusterings done in the order they occured
   */
  list<ClusteringPtr> _clusterings;

  /**
   * The current particle content.
   */
  vector<tClusteringParticlePtr> _currentParticles;

  /**
   * The current guide
   */
  tClusteringGuidePtr _guide;

  /**
   * The hard process object which terminated the clustering
   */
  tCKKWHardProcessPtr _hard;

  /**
   * A stack of clusterings and associated
   * following guides -- ordered from ClusteringSelector
   */
  list<list<pair<ClusteringPtr, tClusteringGuidePtr> > > _allClusterings;

  /**
   * Conversion map of ThePEG::Particle
   */
  map<tPPtr,tClusteringParticlePtr> _conversionMap;

  //@}

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<CascadeReconstructor> initCascadeReconstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CascadeReconstructor & operator=(const CascadeReconstructor &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CascadeReconstructor. */
template <>
struct BaseClassTrait<Herwig::CascadeReconstructor,1> {
  /** Typedef of the first base class of CascadeReconstructor. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CascadeReconstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::CascadeReconstructor>
  : public ClassTraitsBase<Herwig::CascadeReconstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::CascadeReconstructor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * CascadeReconstructor is implemented. It may also include several, space-separated,
   * libraries if the class CascadeReconstructor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "CascadeReconstructor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CascadeReconstructor.tcc"
#endif

#endif /* HERWIG_CascadeReconstructor_H */
