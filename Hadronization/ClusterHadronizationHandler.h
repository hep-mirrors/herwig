// -*- C++ -*-
#ifndef HERWIG_ClusterHadronizationHandler_H
#define HERWIG_ClusterHadronizationHandler_H

#include <ThePEG/Handlers/HadronizationHandler.h>
#include "PartonSplitter.h"
#include "ClusterFinder.h"
#include "ColourReconnector.h"
#include "ClusterFissioner.h"
#include "LightClusterDecayer.h"
#include "ClusterDecayer.h"
#include "ForcedSplitting.h"
#include "Cluster.h"
#include "ClusterHadronizationHandler.fh"

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

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ClusterHadronizationHandler();
  //@}

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

  /**
   * It returns minimum virtuality^2 of partons to use in calculating 
   * distances. It is used both in the Showering and Hadronization.
   */
  inline Energy2 minVirtuality2() const;

  /**
   * It returns the maximum displacement that is allowed for a particle
   * (used to determine the position of a cluster with two components).
   */
  inline Length maxDisplacement() const;

  /**
   * It returns true/false according if the soft underlying model
   * is switched on/off. 
   */
  inline bool isSoftUnderlyingEventON() const;

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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

   /**
    * Initialize this object at the begining of the run phase.
    */
  virtual void doinitrun();
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

  /**
   * The minimum virtuality^2 of partons to use in calculating 
   * distances.
   */
  Energy2 _minVirtuality2;

  /**
   * The maximum displacement that is allowed for a particle
   * (used to determine the position of a cluster with two components).
   */
  Length _maxDisplacement;

  /**
   *  Is the soft underlying event on/off
   */
  bool _softUnderlyingEventMode;

  /**
   * The pointer to the Underlying Event handler. If _softUnderlyingEventMode is true, 
   * this pointer must be set. This is checked in doinit().
   */
  StepHdlPtr _underlyingEventHandler;
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
  static string library() { return ""; }
};

}

#include "ClusterHadronizationHandler.icc"

#endif /* HERWIG_ClusterHadronizationHandler_H */
