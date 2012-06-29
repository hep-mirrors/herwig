// -*- C++ -*-
//
// ClusterHadronizationHandler.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ClusterHadronizationHandler_H
#define HERWIG_ClusterHadronizationHandler_H

#include <ThePEG/Handlers/HadronizationHandler.h>
#include "PartonSplitter.h"
#include "ClusterFinder.h"
#include "ColourReconnector.h"
#include "ClusterFissioner.h"
#include "LightClusterDecayer.h"
#include "ClusterDecayer.h"
#include "ClusterHadronizationHandler.fh"

namespace Herwig {
using namespace ThePEG;


/** \ingroup Hadronization
 *  \class ClusterHadronizationHandler
 *  \brief Class that controls the cluster hadronization algorithm.
 *  \author Philip Stephens  //  cerr << *ch.currentEvent() << '\n';
  cerr << finalHadrons.size() << '\n';

  cerr << "Finished hadronizing \n";

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
 * @see \ref ClusterHadronizationHandlerInterfaces "The interfaces"
 * defined for ClusterHadronizationHandler.
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
		      const Hint & hint);

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
   * The pointer to the Underlying Event handler. 
   */
  StepHdlPtr _underlyingEventHandler;
};


}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

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
  static string className() { return "Herwig::ClusterHadronizationHandler"; }
};

/** @endcond */

}

#include "ClusterHadronizationHandler.icc"

#endif /* HERWIG_ClusterHadronizationHandler_H */
