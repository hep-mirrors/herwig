// -*- C++ -*-
#ifndef HERWIG_ClusterHadronizationHandler_H
#define HERWIG_ClusterHadronizationHandler_H

#include <ThePEG/Handlers/HadronizationHandler.h>
#include "Herwig++/Utilities/GlobalParameters.h"
#include "PartonSplitter.h"
#include "ClusterFinder.h"
#include "ColourReconnector.h"
#include "ClusterFissioner.h"
#include "LightClusterDecayer.h"
#include "ClusterDecayer.h"
#include "ForcedSplitting.h"
#include "Cluster.h"


namespace Herwig {
using namespace ThePEG;


/** \ingroup Hadronization
 *  \class ClusterHadronizationHandler
 *  \brief Class that controls the cluster hadronization algorithm.
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *
 *  This class is the main driver of the cluster hadronization: it is 
 *  responsible for the proper handling of all other specific collaborating
 *  classes PartonSplitter, ClusterFinder, ColourReconnector, ClusterFissioner, 
 *  LightClusterDecayer, ClusterDecayer; 
 *  and for the storing of the produced particles in the Event record.
 * 
 *  Notice that the access to the GlobalParameters class 
 *  instance is provided only to allow non-interfaced and non-persistent classes
 *  (Cluster) to access the global parameters and/or to draw 
 *  random numbers. This is done in the run initialization, doinitrun()
 *  by setting static pointers defined in those non-interfaced and 
 *  non-persistent classes.
 *
 *  @see GlobalParameters
 *  @see PartonSplitter
 *  @see ClusterFinder
 *  @see ColourReconnector
 *  @see ClusterFissioner
 *  @see LightClusterDecayer
 *  @see ClusterDecayer
 *  @see Cluster
 */ 
class ClusterHadronizationHandler: public HadronizationHandler {

public:

  /**
   * The main method which manages the all cluster hadronization.
   *
   * This routine directs "traffic". It determines which function is called
   * and on which particles/clusters. This function also handles the 
   * situation of vetos on the hadronization.
   */
  virtual void handle(EventHandler & ch, const tPVector & tagged,
		      const Hint & hint) throw(Veto, Stop, Exception);

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object to the begining of the run phase.
   */
  inline virtual void doinitrun();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ClusterHadronizationHandler> initClusterHadronizationHandler;

  /**
   * Private and non-existent assignment operator.
   */
  ClusterHadronizationHandler & operator=(const ClusterHadronizationHandler &);

  /**
   * Print the step for debugging.
   */
  void printStep(tStepPtr ptrStep, const string & title);

  /**
   * Print information about the final, complete collections of clusters.
   */
  void debuggingInfo(EventHandler & ch, ClusterVector &);

  /**
   * This is a pointer to a Herwig::GlobalParameters object.
   */
  GlobParamPtr           _globalParameters;

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
   * This is a poniter to the Herwig::ForcedSplitting object
   */
  ForcedSplittingPtr _forcedSplitter;
};


}

namespace ThePEG {

template <>
/**
 * The following template specialization informs ThePEG about the
 * base class of ClusterHadronizationHandler.
 */
struct BaseClassTrait<Herwig::ClusterHadronizationHandler,1> {
  /** Typedef of the base class of ClusterHadronizationHandler. */
  typedef HadronizationHandler NthBase;
};

template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
struct ClassTraits<Herwig::ClusterHadronizationHandler>: 
    public ClassTraitsBase<Herwig::ClusterHadronizationHandler> {
  /** Return the class name.*/
  static string className() { return "Herwig++::ClusterHadronizationHandler"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwHadronization.so"; }
};

}

#include "ClusterHadronizationHandler.icc"

#endif /* HERWIG_ClusterHadronizationHandler_H */
