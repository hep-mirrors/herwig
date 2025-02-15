// -*- C++ -*-
//
// ClusterDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ClusterDecayer_H
#define HERWIG_ClusterDecayer_H

#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/EventRecord/Step.h>
#include "CluHadConfig.h"
#include "HadronSpectrum.h"
#include "ClusterDecayer.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Hadronization
 *  \class ClusterDecayer
 *  \brief This class decays the "normal" clusters
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *
 *  This class decays the "normal" clusters, e.g. ones that are not heavy
 *  enough for fission, and not too light to decay into one hadron.
 *
 *  This class is directs the production of hadrons via 2-body cluster decays.
 *  The selection of the hadron flavours is given by Herwig::HadronSpectrum.
 *
 *  @see HadronSpectrum
 * @see \ref ClusterDecayerInterfaces "The interfaces"
 * defined for ClusterDecayer.
 */
class ClusterDecayer: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  ClusterDecayer();
  //@}

  /** Decays all remaining clusters into hadrons.
   * This routine decays the clusters that are left after
   * Herwig::ClusterFissioner::fission and
   * Herwig::LightClusterDecayer::decay have been called. These are all
   * the "normal" clusters which are not forced into hadrons by
   * the other functions.
   */
  void decay(const ClusterVector & clusters, tPVector & finalhadrons)
   ;

  /**
   * Return the hadron spectrum
   */
  Ptr<HadronSpectrum>::tptr spectrum() const {
    return _hadronSpectrum;
  }

public:

  /**
   * Standard ThePEG function for writing a persistent stream.
   */
  void persistentOutput(PersistentOStream &) const;

  /**
   * Standard ThePEG function for reading from a persistent stream.
   */
  void persistentInput(PersistentIStream &, int);

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

//@}

private:

  /**
   * Private and non-existent assignment operator.
   */
  ClusterDecayer & operator=(const ClusterDecayer &) = delete;

public:

  /** Decays the cluster into two hadrons.
   *
   *  This routine is used to take a given cluster and decay it into
   *  two hadrons which are returned. If one of the constituents is from
   *  the perturbative regime then the direction of the perturbative parton
   *  is remembered and the decay is preferentially in that direction. The
   *  direction of the decay is given by
   *  \f[ \cos \theta = 1 + S \log r_1 \f]
   *  where \f$ S \f$ is a parameter of the model and \f$ r_1 \f$ is a random
   *  number [0,1].
   */
  pair<PPtr,PPtr> decayIntoTwoHadrons(tClusterPtr ptr);
  pair<PPtr,PPtr> decayIntoTwoHadronsDefault(tClusterPtr ptr);
  pair<PPtr,PPtr> decayIntoTwoHadronsNew(tClusterPtr ptr);

private:

  /** Compute the positions of the new hadrons based on the clusters position.
   *
   *  This method calculates the positions of the children hadrons by a
   *  call to ThePEG::RandomGenerator::rndGaussTwoNumbers with width inversely
   *  proportional to the cluster mass, around the parent cluster position.
   */
  void calculatePositions( const Lorentz5Momentum &, const LorentzPoint &,
			   const Lorentz5Momentum &, const Lorentz5Momentum &,
			   LorentzPoint &, LorentzPoint &) const;

  //@{
  /**
   * Whether a cluster decays along the perturbative parton direction.
   */
  bool _clDirLight;
  map<long,bool> _clDirHeavy;
  bool _clDirExotic;

   /**
   * The S parameter from decayIntoTwoHadrons
   */
  double _clSmrLight;
  map<long,double> _clSmrHeavy;
  double _clSmrExotic;
  //@}

  /**
   * Whether or not the hadrons produced should be on-shell
   * or generated used the MassGenerator
   */
  bool _onshell;

  /**
   * Number of tries to generate the masses of the decay products
   */
  unsigned int _masstry;

  /**
   * Use new Kinematics
   */
  int _kinematics;

  /**
   * The hadron spectrum to consider
   */
  Ptr<HadronSpectrum>::ptr _hadronSpectrum;

};


}

#endif /* HERWIG_ClusterDecayer_H */
