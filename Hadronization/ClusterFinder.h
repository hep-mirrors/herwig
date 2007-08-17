// -*- C++ -*-
#ifndef HERWIG_ClusterFinder_H
#define HERWIG_ClusterFinder_H

#include <ThePEG/Interface/Interfaced.h>
#include "CluHadConfig.h"
#include "ClusterFinder.fh"

namespace Herwig {
using namespace ThePEG;

/*! \ingroup Hadronization
 *  \class ClusterFinder
 *  \brief This class forms clusters from the partons produced in the Shower.
 *  \author Philip Stephens
 *  \author Alberto Ribon
 * 
 *  This class scans through the particles in the event and produces a 
 *  collection of clusters, defined as a colour-singlet combinations of 
 *  colour-connected particles. There are no assumptions about the type 
 *  (i.e. quark or diquark) or number of the component particles of the 
 *  cluster; however, most of the time clusters are formed by quark-antiquark 
 *  pairs. In special situations, such as baryon-violating processes in
 *  R-nonconserved Susy, three quarks (or three antiquarks) could form a 
 *  cluster. Because at the moment we don't know how to handle 3-component 
 *  clusters (i.e. how to fission heavy ones, or how to decay clusters), we 
 *  provide also a separate method, reduceToTwoComponents, which 
 *  does the job of redefining these 3-component clusters as "normal" 
 *  2-component ones, simply by randomly considering two (anti-) quarks as a 
 *  (anti-) diquark. Notice that if in the future the method 
 *  reduceToTwoComponents is modified or even eliminated, the 
 *  main method for finding clusters, formClusters, will not need 
 *  any change. 
 *
 * @see \ref ClusterFinderInterfaces "The interfaces"
 * defined for ClusterFinder.
 */
class ClusterFinder: public Interfaced {

public:

  /** 
   * This routine forms the clusters of the event.
   *
   * Form clusters starting from the list of partons given.
   * It also checks if the cluster is a beam cluster, that is if
   * at least one of its components is a beam remnant.
   */
  ClusterVector formClusters(const PVector & partons) 
    throw(Veto, Stop, Exception);

  /**
   * Reduces three component clusters into two components.
   *
   * For the eventual clusters that have three components 
   * (quark, quark, quark) or (antiquark, antiquark, antiquark),
   * it redefines them as "normal" clusters with two components:
   * (quark,diquark) or (antiquark,antidiquark), by a random drawing.
   * This could be eliminated or changed in the future.
   */
  void reduceToTwoComponents(ClusterVector&) 
    throw(Veto, Stop, Exception);

public:

  /**
   * Standard Init function used to initialize the interfaces.
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

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static NoPIOClassDescription<ClusterFinder> initClusterFinder;

  /**
   * Private and non-existent assignment operator.
   */
  ClusterFinder & operator=(const ClusterFinder &);

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

template <>
/**
 * The following template specialization informs ThePEG about the
 * base class of ClusterFinder.
 */
struct BaseClassTrait<Herwig::ClusterFinder,1> {
  /** Typedef of the base class of ClusterFinder. */
  typedef Interfaced NthBase;
};

template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
struct ClassTraits<Herwig::ClusterFinder>: public ClassTraitsBase<Herwig::ClusterFinder> {
  /** Return the class name.*/
  static string className() { return "Herwig::ClusterFinder"; }
};

/** @endcond */

}

#include "ClusterFinder.icc"

#endif /* HERWIG_ClusterFinder_H */
