// -*- C++ -*-
#ifndef HERWIG_PartonicHadronizer_H
#define HERWIG_PartonicHadronizer_H
//
// This is the declaration of the PartonicHadronizer class.
//

#include <ThePEG/Handlers/EventHandler.h>
#include "PartonSplitter.h"
#include "ClusterFinder.h"
#include "ColourReconnector.h"
#include "ClusterFissioner.h"
#include "LightClusterDecayer.h"
#include "ClusterDecayer.h"
#include "Cluster.h"
#include "ThePEG/Interface/Interfaced.h"
#include "PartonicHadronizer.fh"
#include "Cluster.h"  

namespace Herwig {
using namespace ThePEG;

/**
 *  Define some sets
 */
ThePEG_DECLARE_MULTISET(tcPDPtr,cParticleMSet);

/** \ingroup Hadronization
 *  \class PartonicHadronizer
 *  \brief Performs cluster hadronization in partonic bottom and charm meson decays.
 *  \author Peter Richardson
 *
 *  This class is based on the ClusterHadronization handler and is designed
 *  to be called by the HwDecayHandler to perform the hadronization of partonic
 *  bottom and charm meson decays. This is done so that these decays can be vetoed
 *  if the hadronization produces a decay mode which is all ready included as an
 *  inclusive mode to avoid double counting.
 * 
 *  @see PartonSplitter
 *  @see ClusterFinder
 *  @see ColourReconnector
 *  @see ClusterFissioner
 *  @see LightClusterDecayer
 *  @see ClusterDecayer
 *  @see Cluster
 *  @see ClusterHadronizationHandler
 */ 
class PartonicHadronizer: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline PartonicHadronizer();
  //@}


public:

  /**
   * The main method which manages the all cluster hadronization.
   *
   * This routine directs "traffic". It determines which function is called
   * and on which particles/clusters. This function also handles the 
   * situation of vetos on the hadronization.
   * @param parent The parent for the decay.
   * @param pstep The step in which to place the products.
   * @param ch The event handler
   * @param hadrons The hadrons produced in the decay
   * @return Whether or not the hadronization was successful.
   */
  bool hadronize(tPPtr parent,StepPtr pstep,EventHandler & ch,vector<tPPtr> & hadrons);

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<PartonicHadronizer> initPartonicHadronizer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PartonicHadronizer & operator=(const PartonicHadronizer &);

private:

  /** @name Members to handle the vetoing of decays where the cluster model
   * reproduces inclusive modes. 
   */
  //@{

  /**
   * Find the clusters produced in partonic hadron decays
   * @param pstep The step to search
   * @param parent The decaying particles.
   * @param clusters The clusters produced in the decay.
   */
  void findPartonicClusters(Step & pstep,tPPtr parent,vector<tcPPtr> & clusters);
  
  /**
   * Does a cluster only decay to hadrons
   * @param clu The cluster
   * @return Whether or not the cluster decays to hadrons.
   */
  bool hadronicCluster(tPPtr clu);

  /**
   * Check hadrons produce in a partonic hadron decay do not reproduce an inclusive
   * mode.
   * @param parent The decaying particles.
   * @param clusters The clusters produced in the decay.
   * @param hadrons The hadrons produced in the partonic decay.
   * @return Whether or not there are duplicate modes.
   */
  bool duplicateMode(tPPtr parent,vector<tcPPtr> & clusters,vector<tPPtr> & hadrons);

  /**
   *  Remove quarks from cluster splitting from the event record
   * @param outhad Outgoing hadron
   * @param quarks The quarks to be removed
   */
  void removeQuarks(tPPtr outhad,ParticleVector & quarks);
  //@}

private:

  /**
   * This is a pointer to a Herwig::PartonSplitter object.
   */
  PartonSplitterPtr      _partonSplitter;

  /**
   * This is a pointer to a Herwig::ClusterFinder object.
   */
  ClusterFinderPtr       _clusterFinder;

  /**
   * This is a pointer to a Herwig::ColourReconnector object.
   */
  ColourReconnectorPtr   _colourReconnector;

  /**
   * This is a pointer to a Herwig::ClusterFissioner object.
   */
  ClusterFissionerPtr    _clusterFissioner;

  /**
   * This is a pointer to a Herwig::LightClusterDecayer object.
   */
  LightClusterDecayerPtr _lightClusterDecayer;

  /**
   * This is a pointer to a Herwig::ClusterDecayer object.
   */
  ClusterDecayerPtr      _clusterDecayer; 

  /**
   * Switch to control hadrons produced in partonic b and c decays
   */
  bool _exclusive;

  /**
   *  Number of tries for partonic modes
   */
  unsigned int _partontries;

  /**
   * Whether or not to include the intermediates
   */
  bool _inter;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PartonicHadronizer. */
template <>
struct BaseClassTrait<Herwig::PartonicHadronizer,1> {
  /** Typedef of the first base class of PartonicHadronizer. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the PartonicHadronizer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::PartonicHadronizer>
  : public ClassTraitsBase<Herwig::PartonicHadronizer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::PartonicHadronizer"; }
};

/** @endcond */

}

#include "PartonicHadronizer.icc"

#endif /* HERWIG_PartonicHadronizer_H */
