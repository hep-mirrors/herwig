// -*- C++ -*-
//
// PartonicDecayerBase.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PartonicDecayerBase_H
#define HERWIG_PartonicDecayerBase_H
//
// This is the declaration of the PartonicDecayerBase class.
//

#include "Herwig/Decay/HwDecayerBase.h"
#include "Herwig/Hadronization/PartonSplitter.h"
#include "Herwig/Hadronization/ClusterFinder.h"
#include "Herwig/Hadronization/ClusterFissioner.h"
#include "Herwig/Hadronization/LightClusterDecayer.h"
#include "Herwig/Hadronization/ClusterDecayer.h"
#include "Herwig/Hadronization/Cluster.h"

namespace Herwig {

using namespace ThePEG;

/**
 *  Define some sets
 */
ThePEG_DECLARE_MULTISET(tcPDPtr,cParticleMSet);

/**
 * Here is the documentation of the PartonicDecayerBase class.
 *
 * @see \ref PartonicDecayerBaseInterfaces "The interfaces"
 * defined for PartonicDecayerBase.
 */
class PartonicDecayerBase: public HwDecayerBase {

public:

  /**
   * The default constructor.
   */
  PartonicDecayerBase();

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay for a particular mode
   * @param parent The decaying particle
   * @param children The decay products
   * @return true If this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const = 0;
  
  /**
   *  Perform the decay of the particle to the specified decay products
   * @param parent The decaying particle
   * @param children The decay products
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const tPDVector & children) const = 0;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;
  //@}

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

  /**
   * Check hadrons produced in a partonic hadron decay do not reproduce an inclusive
   * mode.
   * @param parent The decaying particles.
   * @param hadrons The hadrons produced in the partonic decay.
   * @return Whether or not there are duplicate modes.
   */
  bool duplicateMode(const Particle & parent,
		     const vector<tPPtr> & hadrons) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PartonicDecayerBase & operator=(const PartonicDecayerBase &) = delete;

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

#endif /* HERWIG_PartonicDecayerBase_H */
