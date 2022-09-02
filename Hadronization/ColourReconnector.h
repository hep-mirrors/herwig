// -*- C++ -*-
//
// ColourReconnector.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
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
   * @brief     calculates the "euclidean" distance of two quarks in the 
   * 			rapdidity-phi plane
   * @arguments p, q the two quarks
   * @return    the dimensionless distance:
   * 			\deltaR_12=sqrt(\deltaY_12^2 + \delta\phi_12^2)
   */
  double _displacement(tcPPtr p, tcPPtr q) const;
  
   
  /**
   * @brief     calculates the "euclidean" distance of a the 3 (anti)quarks
   *			(anti)baryonic cluster in the rapdidity-phi plane
   * @arguments p, q the two quarks
   * @return    the dimensionless distance:
   * 			if Junction is enabled the difference of all 3 quarks
   * 			with respect to their mean point is calculated:
   * 				<Y>    = sum_i Y_i/3
   * 				<\phi> = sum_i \phi_i/3
   * 				\deltaR     = sum_i sqrt( (Y_i - <Y>)^2 + (\phi_i - <phi>)^2)
   *
   * 			if Junction is disabled the difference of all distinct
   * 			pairing of the 3 quarks is computed:
   * 				\deltaR_ij  = sqrt(\deltaY_ij^2 + \delta\phi_ij^2)
   * 				\deltaR_tot = sum_{i,j<i} \deltaR_ij
   * 			
   * 			NOTE:   switch the Junction option will necessarily
   * 					need to be combined with a change (or re-tune)
   * 					of _mesonToBaryonFactor otherwise no well 
   * 					description is to be expected
   * 			TODO:	maybe add a different p-norm option to get
   * 					more phenomenology
   */
  double _displacementBaryonic(tcPPtr p1, tcPPtr p2, tcPPtr p3) const;

  
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
  
  
  /**
   * @brief     BaryonicMesonic colour reconnection model
   * @arguments cv cluster vector
   * BaryonicMesonic Colour Reconnection model with reconnections from mesonic clusters
   * to baryonic cluster and the contrary. Allows also reconnection between mesonic
   * and baryonic Clusters 
   */
  void _doRecoBaryonicMesonic(ClusterVector & cv) const; 
  
  


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
   * @brief     Finds best CR option for the BaryonicMesonic CR model
   * @arguments cls vector of selected clusters, baryonic is the number of baryonic
   * 			clusters in selection, kind_of_reco is output string denoting best
   * 			CR option and reco_info is output storage for information on which
   * 			(anti-)quarks to exchange and in which way.
   * BaryonicMesonic Colour Reconnection model with reconnections from mesonic clusters
   * to baryonic cluster and the contrary. Allows also reconnection between mesonic
   * and baryonic Clusters 
   */
  void _findbestreconnectionoption(std::vector<CluVecIt> &cls,
		  									const unsigned &baryonic,
											string &kind_of_reco,
											int (&reco_info)[3]) const;

  /**
   * @brief     Reconnects the constituents of the given clusters to the (only)
   *            other possible cluster combination.
   * @return    pair of pointers to the two new clusters
   * Used for Plain and Baryonic Colour Reconnection models
   */
  pair <ClusterPtr,ClusterPtr> _reconnect(ClusterPtr &c1, ClusterPtr &c2) const;
  
  /**
   * @brief     Reconnects (2B->2B) the constituents of the given clusters to 
    			another possible cluster combination whose information is given 
    			in s1 and s2.
   * @arguments c1 and c2 are baryonic clusters and s1 and s2 are the respective
   				indices of the constituents of c1 and c2 respectively 
   * @return    pair of pointers to the two new clusters
   * Used only in BaryonicMesonic algorithm and will exchange constituent s1 of
   * c1 with constituent s2 of c2
   */
  pair <ClusterPtr,ClusterPtr> _reconnect2Bto2B(ClusterPtr &c1, ClusterPtr &c2, const int s1, const int s2) const;
  
  /**
   * @brief     Reconnects (B,Bbar->3M) the constituents of the given clusters to 
   				another	possible cluster combination whose information is given in 
				s0, s1 and s2.
   * @arguments c1 and c2 are baryonic clusters. s0, s1 and s2 are the respective
   				indices which determine the CR 
   * @return    tuple of pointers to the 3 new mesonic clusters
   * Used only in BaryonicMesonic algorithm and will form 3 mesonic clusters according
   * to the indices s0, s1 and s2. The i-th constituent of c1 is connected to the si-th 
   * constituent of c2
   */
  std::tuple  <ClusterPtr, ClusterPtr, ClusterPtr> _reconnectBbarBto3M(ClusterPtr &c1, ClusterPtr &c2, const int s0, const int s1, const int s2 ) const;
  
  /**
   * @brief     Reconnects (3M->3M) the constituents of the given clusters to 
   				another	possible cluster combination whose information is given in 
				sinfos
   * @arguments c1, c2 and c3 are mesonic clusters. infos[3] are the respective
   				indices which determine the CR 
   * @return    tuple of pointers to the 3 CR'd mesonic clusters
   * Used only in BaryonicMesonic algorithm and will reconnect 3 mesonic clusters according
   * to the infos, which determine the CR. The coloured quark of the i-th cluster forms 
   * a new cluster with the anticoloured quark of the info[i]-th cluster
   */
  std::tuple  <ClusterPtr, ClusterPtr, ClusterPtr> _reconnect3Mto3M(ClusterPtr &c1, ClusterPtr &c2, ClusterPtr &c3, const int infos[3] ) const;
  
   
    /**
   * Reconnection method for a Baryonic and a Mesonic Cluster to a Baryonic and a Mesonic Cluster
   * s0 is the Number of the (Anti)Patrticle of the Baryonic Cluster , which should be swapped with the Anti(Particle) of the Mesonic Cluster
   */ 
  /**
   * @brief     Reconnects the constituents of the given clusters to 
   				another	possible cluster combination whose information 
				is given in s0.
   * @arguments c1 and c2 are one baryonic and one mesonic cluster respectively 
   				and s0 is the respective index of the constituents of the baryonic
				cluster which is to be exchangeed with the mesonic cluster.
   * @return    pair of pointers to the two new clusters
   * Used only in BaryonicMesonic algorithm and will exchange constituent s0 of
   * the baryonic cluster with the (anti)-quark of the mesonic cluster
   */
  pair  <ClusterPtr, ClusterPtr> _reconnectMBtoMB(ClusterPtr &c1, ClusterPtr &c2, const int s0) const;
  
  
  
  
  /**
   * Methods for the BaryonicMesonic CR algorithm
   * Find the best reconnection option for the respective cluster-combination 
   * 
   */

  /**
   * @brief 	veto for far apart clusters
   * @arguments	expects at most 3 CluVecIt in clu vector
   * @return	returns true if clusters are more distant than _maxDistance 
   * 			in space
   * TODO: problematic maybe add option to turn off
   */
  bool _clustersFarApart( const std::vector<CluVecIt> & clu ) const;
	
  /**
   * @brief     Does reconnect 2 beam clusters for BaryonicMesonic CR model
   * 			if option CR2BeamClusters is enabled
   * @arguments cv cluster vector
   */
  void _doReco2BeamClusters(ClusterVector & cv) const;


  /**
   * @brief     finds the best reconnection option and stores it in bswap1
   * 			and bswap2 (2B->2B colour reconnection)
   * @arguments c1 and c2 cluster pointers and kind_of_reco will output
   * 			the best found reconnection option for c1 and c2
   */
  void _2Bto2BreconnectionFinder(ClusterPtr &c1, ClusterPtr &c2, int &bswap1, int &bswap2, double mindisplsum, string &kind_of_reco ) const;

  /**
   * @brief     finds the best reconnection option and stores it in mswap0
   * 			mswap1 and mswap2 (BbarB->3M colour reconnection)
   * @arguments c1 and c2 cluster pointers and kind_of_reco will output
   * 			the best found reconnection option for c1 and c2
   */
  void _BbarBto3MreconnectionFinder(ClusterPtr &c1, ClusterPtr &c2, int &mswap0, int &mswap1, int &mswap2,
															double min_displ_sum, string & kind_of_reco) const; 

  /**
   * @brief     finds the best reconnection option and stores it in swap
   * 			(BM->BM colour reconnection)
   * @arguments c1 and c2 cluster pointers and kind_of_reco will output
   * 			the best found reconnection option for c1 and c2
   */
  void _BMtoBMreconnectionfinder(ClusterPtr &c1, ClusterPtr &c2, int &swap, 
															double min_displ_sum, string &kind_of_reco) const;

  /**
   * @brief     finds the best reconnection option and stores it in swap0,
   * 			swap1 and swap2	(3M->{3M,BbarB} colour reconnection)
   * @arguments c1 and c2 cluster pointers and kind_of_reco will output
   * 			the best found reconnection option for c1 and c2
   */
  void _3MtoXreconnectionfinder(std::vector<CluVecIt> &cv, int &swap0, int &swap1,
																int &swap2, double min_displ_sum, string &kind_of_reco) const;

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

private :

  /** Data Members */
  //@
  /**
   * Specifies the colour reconnection algorithm to be used.
   */
  int _algorithm = 0;

  /**
   * Do we do colour reconnections?
   */
  int _clreco = 0;

  /**
   * Do we want debug informations? 
   */
  int _debug = 0;
  
  /**
   * @brief     Junction-like model for BaryonicMesonic model
   * Do we want to use the junction-like model for
   * computing the displacements of BaryonicMesonic model?
   * otherwise pairwise distances are used.
   * If _junctionMBCR is activated the displacements are computed in the
   * rapidity-Phi plane by difference to the average rapidity and phi:
   * DeltaR_i^2 = (rapidity_i - meanRap)^2 + (phi_i - meanPhi)^2 
   * DeltaR = sum_i DeltaR_i 
   * if _junctionMBCR=0 the displacements are computed:
   * DeltaR_ij^2 = (rapidity_i - rapidity_j)^2 + (phi_i - phi_j)^2 
   * DeltaR = sum_i,j<i DeltaR_ij
   */
  int _junctionMBCR = 1;

  /**
   * Statistical Reco: 
   * Factor used to determine the initial temperature according to
   * InitialTemperature = _initTemp * median {energy changes in a few random
   * rearrangements}
   */
  double _initTemp = 0.01;

  /**
   * Statistical Reco: 
   * The annealing factor is the ratio of two successive temperature steps:
   * T_n = _annealingFactor * T_(n-1)
   */
  double _annealingFactor = 0.21;

  /**
   * Statistical Reco: 
   * Number of temperature steps in the statistical annealing algorithm
   */
  unsigned int _annealingSteps = 10;

  /**
   * Statistical Reco:
   * The number of tries per temperature steps is the number of clusters times
   * this factor.
   */
  double _triesPerStepFactor = 0.66;

  /**
   * Probability that a found reconnection possibility is actually accepted.
   * used in Plain & Baryonic CR
   */
  double _preco = 0.5;

  /**
   * Probability that a found reconnection possibility is actually accepted.
   * used in Baryonic CR
   */
  double _precoBaryonic = 0.5;

  /**
   * Probability that a found reconnection possibility is actually accepted.
   * For reconnecting 3M -> 3M'
   * used in BaryonicMesonic
   * NOTE: if 0 this type of reconnection is not even tried
   */
  double _preco3M_3M = 0.5;

  /**
   * Probability that a found reconnection possibility is actually accepted.
   * For reconnecting 3M -> B,Bbar
   * used in BaryonicMesonic
   */
  double _preco3M_BBbar = 0.5;

  /**
   * Probability that a found reconnection possibility is actually accepted.
   * For reconnecting Bbar,B -> 3M
   * used in BaryonicMesonic
   */
  double _precoBbarB_3M = 0.5;
  
  /**
   * Probability that a found reconnection possibility is actually accepted.
   * For reconnecting 2B -> 2B' or 2Bbar -> 2Bbar' 
   * used in BaryonicMesonic
   * NOTE: if 0 this type of reconnection is not even tried
   */
  double _preco2B_2B = 0.5;
  
  /**
   * Probability that a found reconnection possibility is actually accepted.
   * For reconnecting M,B -> M',B' or M,Bbar -> M',Bbar'
   * used in BaryonicMesonic
   * NOTE: if 0 this type of reconnection is not even tried
   */
  double _precoMB_MB = 0.5;
  
  /**
   * For the BaryonicMesonic algorithm
   * How many times do suggest cluster for reconnection?
   * n(reconnectionstries) = _stepFactor * n(clusters)*n(clusters);
   */ 
  double _stepFactor = 5.0;

  /**
   * Factor for comparing mesonic clusters to baryonic clusters
   */
  double _mesonToBaryonFactor = 2.0;

  /**
   *  Maximium distance for reconnections
   *  TODO: remove if issues with anticausality are discussed and resolved
   */
  Length _maxDistance = femtometer;

  /**
   * @return	true, if the two partons are splitting products of the same
   * 		gluon
   */
  bool _isColour8(tcPPtr p, tcPPtr q) const;
  
  /**s
   *  Option for handling octets
   */
  unsigned int _octetOption = 0;

  /**
   *  Option for colour reconnecting 2 Beam Clusters if no others are present
   */
  int _cr2BeamClusters = 0;

  /**
   *  Option for performing Plain colour reconnection before the Statistical,
   *  Baryonic or BaryonicMesonic algorithm is performed
   */
  int _prePlainCR = 0;

  /**
   *  Option for colour reconnecting Clusters only if their vertex 3-distance
   *  is less than _maxDistance
   */
  int _localCR = 0;

  /**
   *  Option for colour reconnecting Clusters only if their spacetime difference
   *  is bigger than 0
   */
  int _causalCR = 0;

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   *   Option for doubly heavy baryons
   */
  bool _doublyHeavyBaryon = true;
  //@}
};


}

#endif /* HERWIG_ColourReconnector_H */

