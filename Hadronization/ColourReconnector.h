// -*- C++ -*-
//
// ColourReconnector.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ColourReconnector_H
#define HERWIG_ColourReconnector_H

#include <ThePEG/Interface/Interfaced.h>
#include "CluHadConfig.h"
#include "ColourReconnector.fh"

namespace Herwig {


using namespace ThePEG;

/** \ingroup Hadronization
 *  \class ColourReconnector
 *  \brief Class for changing colour reconnections of partons.
 *  \author Alberto Ribon, Christian Roehr
 * 
 *  This class does the nonperturbative colour rearrangement, after the 
 *  nonperturbative gluon splitting and the "normal" cluster formation. 
 *  It uses the list of particles in the event record, and the collections of
 *  "usual" clusters which is passed to the main method. If the colour 
 *  reconnection is actually accepted, then the previous collections of "usual"
 *  clusters is first deleted and then the new one is created.
 *
 * * @see \ref ColourReconnectorInterfaces "The interfaces"
 * defined for ColourReconnector.
 */
class ColourReconnector: public Interfaced {

public:

  /**
   * Does the colour rearrangement, starting out from the list of particles in
   * the event record and the collection of "usual" clusters passed as
   * arguments. If the actual rearrangement is accepted, the initial collection of
   * clusters is overridden by the old ones.
   */
  void rearrange(ClusterVector & clusters);

  using CluVecIt = ClusterVector::iterator;

private:

  /** PRIVATE MEMBER FUNCTIONS */

  /**
   * @brief     Calculates the sum of the squared cluster masses.
   * @arguments q, aq vectors containing the quarks and antiquarks respectively
   * @return    Sum of cluster squared masses M^2_{q[i],aq[i]}.
   */
  Energy2 _clusterMassSum(const PVector & q, const PVector & aq) const;

  
  /**
   * @brief  Examines whether the cluster vector (under the given permutation of
   *         the antiquarks) contains colour-octet clusters
   * @param  cv Cluster vector
   * @param  P  Permutation, a vector of permutated indices from 0 to
   *            cv.size()-1
   */
  bool _containsColour8(const ClusterVector & cv, const vector<size_t> & P) const;

  /**
   * @brief     A Metropolis-type algorithm which finds a local minimum in the
   *            total sum of cluster masses
   * @arguments cv cluster vector
   */
  void _doRecoStatistical(ClusterVector & cv) const;

  /**
   * @brief     Plain colour reconnection as used in Herwig 2.5.0
   * @arguments cv cluster vector
   */
  void _doRecoPlain(ClusterVector & cv) const;


  /**
   * Baryonic Colour Reconnection model
   */
  void _doRecoBaryonic(ClusterVector & cv) const;


  void _makeBaryonicClusters(ClusterPtr &c1, ClusterPtr &c2, ClusterPtr &c3,
			     ClusterPtr &newcluster1, ClusterPtr &newcluster2) const;
  

  /**
   * @brief     Finds the cluster in cv which, if reconnected with the given
   *            cluster cl, would result in the smallest sum of cluster masses.
   *            If no reconnection partner can be found, a pointer to the
   *            original Cluster cl is returned.
   * @arguments cv cluster vector
   *            cl cluster iterator (must be from cv) which wants to have a reconnection partner
   * @return    iterator to the found cluster, or the original cluster pointer if
   *            no mass-reducing combination can be found
   */


  CluVecIt _findRecoPartner(CluVecIt cl, ClusterVector & cv) const;

  CluVecIt _findPartnerRapidity(CluVecIt cl, ClusterVector & cv) const;

  CluVecIt _findPartnerBaryonic(CluVecIt cl, ClusterVector & cv, 
                                               bool & tetraCand, 
                                               const ClusterVector& a, 
                                               CluVecIt  &baryonic1,
                                               CluVecIt  &baryonic2 ) const;




  /**
   * @brief     Reconnects the constituents of the given clusters to the (only)
   *            other possible cluster combination.
   * @return    pair of pointers to the two new clusters
   */
  pair <ClusterPtr,ClusterPtr> _reconnect(ClusterPtr &c1, ClusterPtr &c2) const;
  
  /**
   * Reconnection method for baryonic reconenction model 
   */ 
  pair <ClusterPtr,ClusterPtr> _reconnectBaryonic(ClusterPtr &c1, ClusterPtr &c2) const;


  /**
   * @brief     At random, swap two antiquarks, if not excluded by the
   *            constraint that there must not be any colour-octet clusters.
   * @arguments q, aq vectors containing the quarks and antiquarks respectively
   *            maxtries  maximal number of tries to find a non-colour-octet
   *                      reconfiguration
   * @return    Pair of ints indicating the indices of the antiquarks to be
   *            swapped. Returns (-1,-1) if no valid reconfiguration could be
   *            found after maxtries trials
   */
  pair <int,int>
    _shuffle(const PVector & q, const PVector & aq, unsigned maxtries = 10) const;


  /** DATA MEMBERS */

  /**
   * Specifies the colour reconnection algorithm to be used.
   */
  int _algorithm = 0;

  /**
   * The annealing factor is the ratio of two successive temperature steps:
   * T_n = _annealingFactor * T_(n-1)
   */
  double _annealingFactor = 0.9;

  /**
   * Number of temperature steps in the statistical annealing algorithm
   */
  unsigned int _annealingSteps = 50;

  /**
   * Do we do colour reconnections?
   */
  int _clreco = 0;

  /**
   * Factor used to determine the initial temperature according to
   * InitialTemperature = _initTemp * median {energy changes in a few random
   * rearrangements}
   */
  double _initTemp = 0.1;

  /**
   * Probability that a found reconnection possibility is actually accepted.
   */
  double _preco = 0.5;


  double _precoBaryonic = 0.5;
  
  /**
   * The number of tries per temperature steps is the number of clusters times
   * this factor.
   */
  double _triesPerStepFactor = 5.0;
  /**
   * maximum allowed distance in the eta phi space for reconnection to occur
   */

  /**
   *  Maximium distance for reconnections
   */
  Length _maxDistance = picometer;

  /**
   * @return	true, if the two partons are splitting products of the same
   * 		gluon
   */
  bool _isColour8(tcPPtr p, tcPPtr q) const;
  
  /**
   *  Option for handling octets
   */
  unsigned int _octetOption = 0;

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
  ColourReconnector & operator=(const ColourReconnector &) = delete;


};


}

#endif /* HERWIG_ColourReconnector_H */

