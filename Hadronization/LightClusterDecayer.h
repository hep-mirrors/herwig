// -*- C++ -*-
//
// LightClusterDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_LightClusterDecayer_H
#define HERWIG_LightClusterDecayer_H

#include <ThePEG/Interface/Interfaced.h>
#include "CluHadConfig.h"
#include "HadronSelector.h"
#include "LightClusterDecayer.fh"


namespace Herwig {


using namespace ThePEG;

/** \ingroup Hadronization
 *  \class LightClusterDecayer
 *  \brief This class performs the decay of light clusters into a single hadron.
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *
 *  This is the class that performs the decay of light clusters into 
 *  only one hadron. The major difficulty is that a kinematical reshuffling
 *  is necessary, between the cluster under consideration and its
 *  "neighbouring" clusters, to conserve energy-momentum in one-body decay.
 *  Notice that, differently from what happens in Fortran Herwig,
 *  light (that is below the threshold for the production of the lightest 
 *  pair of hadrons with the proper flavours) fission products, produced 
 *  by the fission of heavy clusters in ClusterFissioner 
 *  have been already "decayed" into single hadron (the lightest one 
 *  with proper flavour) by the same latter class, without require 
 *  any reshuffling. Therefore the light clusters that are treated in 
 *  this LightClusterDecayer class are produced directly 
 *  (originally) by the ClusterFinder. 
 *	
 *  Notice:
 *  - The choice of the candidate cluster with whom to reshuffle momentum 
 *    is based on the minimal space-time distance from the light cluster.
 *  - An alternate choice of what is considered a "neighbour" could be
 *    implemented but was not considered for Herwig.
 *
 *  @see HadronSelector 
 * @see \ref LightClusterDecayerInterfaces "The interfaces"
 * defined for LightClusterDecayer.
 */ 
class LightClusterDecayer: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  LightClusterDecayer() : _limBottom(), _limCharm(), _limExotic() 
  {} 
  //@}

  /**
   * This method does the decay of light hadron in one hadron.
   *
   * This method requires a kinematical reshuffling for energy-momentum 
   * conservation. This is done explicitly by the (private) method 
   * reshuffling().
   */
  bool decay(ClusterVector & clusters, tPVector & finalhadrons);

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

private:

  /**
   * Private and non-existent assignment operator.
   */
  LightClusterDecayer & operator=(const LightClusterDecayer &);

  /**
   * This (private) method, called by decay(), takes care of the kinematical
   * reshuffling necessary for energy-momentum conservation.
   */
  bool reshuffling( const tcPDPtr, tClusterPtr, tClusterPtr,
		    tClusterVector &, tPVector & finalhadrons) 
   ; 
  
  /**
   *  This (private) method, called by decay(), performs reshuffling in the 
   * special case of a semileptonic partonic b/c decay 
   * @param hadron The hadron to be produced
   * @param cluster The cluster to be reshuffled
   * @param finalhadrons The vector of outgoing hadrons
   */
  bool partonicReshuffle(const tcPDPtr hadron,const PPtr cluster,
			 tPVector & finalhadrons);

  /**
   * A pointer to a Herwig::HadronSelector object used for producing hadrons.
   */
  Ptr<HadronSelector>::pointer _hadronSelector;

  /**
   * @name A parameter used for determining when clusters are too light.
   *
   * This parameter is used for setting the lower threshold, \f$ t \f$ as
   * \f[ t' = t(1 + r B^1_{\rm lim}) \f]
   * where \f$ r \f$ is a random number [0,1].
   */
  //@{
  double _limBottom;
  double _limCharm;
  double _limExotic;
  //@}
};

}

#endif /* HERWIG_LightClusterDecayer_H */
