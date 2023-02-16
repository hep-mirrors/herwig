// -*- C++ -*-
//
// ClusterHadronizationHandler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/Utilities/Reshuffler.h"
#include "GluonMassGenerator.h"

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
class ClusterHadronizationHandler:
    public HadronizationHandler, public Reshuffler {

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
  Energy2 minVirtuality2() const
  { return _minVirtuality2; }

  /**
   * It returns the maximum displacement that is allowed for a particle
   * (used to determine the position of a cluster with two components).
   */
  Length maxDisplacement() const
  { return _maxDisplacement; }

  /**
   * It returns true/false according if the soft underlying model
   * is switched on/off. 
   */
  bool isSoftUnderlyingEventON() const
  { return _underlyingEventHandler; }
  
  /**
   *  pointer to "this", the current HadronizationHandler.
   */
  static const ClusterHadronizationHandler * currentHandler() {
    if(!currentHandler_){
      cerr<< " \nCreating new ClusterHadronizationHandler without input from infiles.";
      cerr<< " \nWhen using for example the string model ";
      cerr<< " hadronic decays are still treated by the Cluster model\n";
      currentHandler_=new ClusterHadronizationHandler();;
    }
    return currentHandler_;
  }

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
  virtual void doinit();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

private:

  /**
   * Private and non-existent assignment operator.
   */
  ClusterHadronizationHandler & operator=(const ClusterHadronizationHandler &) = delete;

  /**
   * This is a pointer to a Herwig::PartonSplitter object.
   */
  map<int,PartonSplitterPtr>      _partonSplitters;

  /**
   * This is a pointer to a Herwig::ClusterFinder object.
   */
  map<int,ClusterFinderPtr>       _clusterFinders;

  /**
   * This is a pointer to a Herwig::ColourReconnector object.
   */
  map<int,ColourReconnectorPtr>   _colourReconnectors;

  /**
   * This is a pointer to a Herwig::ClusterFissioner object.
   */
  map<int,ClusterFissionerPtr>    _clusterFissioners;

  /**
   * This is a pointer to a Herwig::LightClusterDecayer object.
   */
  map<int,LightClusterDecayerPtr> _lightClusterDecayers;
  
  /**
   * This is a pointer to a Herwig::ClusterDecayer object.
   */
  map<int,ClusterDecayerPtr>      _clusterDecayers; 

  /**
   * A pointer to a gluon mass generator for the reshuffling
   */
  map<int,Ptr<GluonMassGenerator>::ptr> _gluonMassGenerators;

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
   * A pointer to a gluon mass generator for the reshuffling
   */
  Ptr<GluonMassGenerator>::ptr _gluonMassGenerator;

  /**
   * Perform reshuffling to constituent masses (0 = no reshuffling, 1 = global reshuffling, 2 = colour connected reshuffling).
   */
  int reshuffle_ = 0;

  /**
   * The minimum virtuality^2 of partons to use in calculating 
   * distances.
   */
  Energy2 _minVirtuality2 = 0.1_GeV2;

  /**
   * The maximum displacement that is allowed for a particle
   * (used to determine the position of a cluster with two components).
   */
  Length _maxDisplacement = 1.0e-10_mm;

  /**
   * The pointer to the Underlying Event handler. 
   */
  StepHdlPtr _underlyingEventHandler;

  /**
   *  How to handle baryon-number clusters
   */
  bool _reduceToTwoComponents = true;

  /**
   * Tag the constituents of the clusters as their parents
   */
  void _setChildren(const ClusterVector & clusters) const;

  /**
   * Split the list of partons into colour connected sub-lists before reshuffling
   */
  void splitIntoColourSinglets(PVector thelist,
			       vector<PVector>& reshufflelists,
			       int interaction);

  /**
   * Use the currently set list of handlers to hadronize the given interaction
   */
  string _setHandlersForInteraction(string); 
  
  /**
   *  pointer to "this", the current HadronizationHandler.
   */
  static ClusterHadronizationHandler * currentHandler_;

};


}

#endif /* HERWIG_ClusterHadronizationHandler_H */
