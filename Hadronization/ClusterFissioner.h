// -*- C++ -*-
//
// ClusterFissioner.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ClusterFissioner_H
#define HERWIG_ClusterFissioner_H

#include <ThePEG/Interface/Interfaced.h>
#include "CluHadConfig.h"
#include "HadronSelector.h"
#include "ClusterFissioner.fh"

namespace Herwig {
using namespace ThePEG;

  //class Cluster;          // forward declaration

/** \ingroup Hadronization
 *  \class ClusterFissioner
 *  \brief This class handles clusters which are too heavy.
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *  \author Stefan Gieseke
 *
 *  This class does the job of chopping up either heavy clusters or beam 
 *  clusters in two lighter ones. The procedure is repeated recursively until 
 *  all of the cluster children have masses below some threshold values.
 *
 *  For the beam remnant clusters, at the moment what is done is the following.
 *  In the case that the soft underlying event is switched on, the 
 *  beam remnant clusters are tagged as not available,
 *  therefore they will not be treated at all during the hadronization. 
 *  In the case instead that the soft underlying event is switched off,
 *  then the beam remnant clusters are treated exactly as "normal" clusters,
 *  with the only exception of the mass spectrum used to generate the
 *  cluster children masses. For non-beam clusters, the masses of the cluster
 *  children are draw from a power-like mass distribution; for beam clusters,
 *  according to the value of the flag _IOpRem, either both 
 *  children masses are draw from a fast-decreasing exponential mass 
 *  distribution (case _IOpRem == 0, or, indendently by 
 *  _IOpRem, in the special case that the beam cluster contains two 
 *  beam remnants), or one mass from the exponential distribution (corresponding
 *  of the cluster child with the beam remnant) and the other with the usual 
 *  power-like distribution (case _IOpRem == 1, which is the 
 *  default one, as in Herwig 6.3). 
 *
 *  The reason behind the use of a fast-decreasing exponential distribution 
 *  is that to avoid a large transverse energy from the many sequential
 *  fissions that would otherwise occur due to the typical large cluster 
 *  mass of beam clusters. Using instead an exponential distribution 
 *  the masses of the two cluster children will be very small (order of 
 *  GeV).
 *
 *  The rationale behind the implementation of the splitting of clusters
 *  has been to preserve *all* of the information about such splitting 
 *  process. More explicitly a ThePEG::Step class is passed in and the
 *  new clusters are added to the step as the decay products of the
 *  heavy cluster. This approach has the twofold 
 *  advantage to provide all of the information that could be needed 
 *  (expecially in future developments), without any information loss, 
 *  and furthermore it allows a better debugging. 
 *
 *  @see HadronSelector 
 * @see \ref ClusterFissionerInterfaces "The interfaces"
 * defined for ClusterFissioner.
 */ 
class ClusterFissioner: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
   ClusterFissioner();
  //@}

  /** Splits the clusters which are too heavy.
   *
   * Split either heavy clusters or beam clusters recursively until all 
   * children have mass below some threshold. Heavy clusters are those that
   * satisfy the condition 
   * \f[ M^P > C^P + S^P \f]
   * where \f$ M \f$ is the clusters mass, \f$ P \f$ is the parameter
   * ClPow, \f$ C \f$ is the parameter ClMax and \f$ S \f$ is the 
   * sum of the clusters constituent partons.
   * For beam clusters, they are split only if the soft underlying event
   * is switched off, otherwise these clusters will be tagged as unavailable
   * and they will not be treated by the hadronization altogether. 
   * In the case beam clusters will be split, the procedure is exactly
   * the same as for normal non-beam clusters, with the only exception
   * of the mass spectrum from which to draw the masses of the two 
   * cluster children (see method drawChildrenMasses for details).
   */
  tPVector fission(ClusterVector & clusters, bool softUEisOn);

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
  ClusterFissioner & operator=(const ClusterFissioner &);

  /** 
   * This method directs the splitting of the heavy clusters
   *
   * This method does the splitting of the clusters and all of its cluster 
   * children, if heavy. All of these new children clusters are added to the
   * collection of clusters. The method works as follows.
   * Initially the vector contains just the stack of input pointers to the
   * clusters to be split. Then it will be filled recursively by all
   * of the cluster's children that are heavy enough to require
   * to be split. In each loop, the last element of the vector is 
   * considered (only once because it is then removed from the vector).
   *
   * \todo is the following still true?
   * For normal, non-beam clusters, a power-like mass distribution
   * is used, whereas for beam clusters a fast-decreasing exponential mass 
   * distribution is used instead. This avoids many iterative splitting which 
   * could produce an unphysical large transverse energy from a supposed 
   * soft beam remnant process.
   */
  void cut(stack<ClusterPtr> &,
	   ClusterVector&, tPVector & finalhadrons, bool softUEisOn);

public:

  /**
   * Definition for easy passing of two particles.
   */
  typedef pair<PPtr,PPtr> PPair;

  /**
   * Definition for use in the cut function.
   */
  typedef pair<PPair,PPair> cutType;

  /** 
   * Splits the input cluster.
   *
   * Split the input cluster (which can be either an heavy non-beam
   * cluster or a beam cluster). The result is two pairs of particles. The
   * first element of each pair is new cluster/hadron, while the second
   * element of each pair is the particle drawn from the vacuum to create
   * the new cluster/hadron.
   * Notice that this method treats also beam clusters by using a different
   * mass spectrum used to generate the cluster child masses (see method
   * drawChildMass).
   */
  //@{
  /**
   *  Split two-component cluster
   */
  virtual cutType cutTwo(ClusterPtr &, tPVector & finalhadrons, bool softUEisOn);

  /**
   *  Split three-component cluster
   */
  virtual cutType cutThree(ClusterPtr &, tPVector & finalhadrons, bool softUEisOn);
  //@}
public:

  /**
   * Produces a hadron and returns the flavour drawn from the vacuum.
   *
   * This routine produces a new hadron. It
   * also sets the momentum and vertex to the values given.
   */
  PPair produceHadron(tcPDPtr hadron, tPPtr newPtr, const Lorentz5Momentum &a,
		      const LorentzPoint &b) const;
protected:

  /**
   * Produces a cluster from the flavours passed in.
   *
   * This routine produces a new cluster with the flavours given by ptrQ and newPtr.
   * The new 5 momentum is a and the parent momentum are c and d. C is for the
   * ptrQ and d is for the new particle newPtr. rem specifies whether the existing
   * particle is a beam remnant or not.
   */
  PPair produceCluster(tPPtr ptrQ, tPPtr newPtr, const Lorentz5Momentum &a,
		       const LorentzPoint &b, const Lorentz5Momentum &c,
		       const Lorentz5Momentum &d, const bool rem,
		       tPPtr spect=tPPtr(), bool remSpect=false) const;

  /**
   * Returns the new quark-antiquark pair
   * needed for fission of a heavy cluster. Equal probabilities
   * are assumed for producing  u, d, or s pairs. 
   */
  void drawNewFlavour(PPtr& newPtrPos,PPtr& newPtrNeg) const;
   
  /**
   * Produces the mass of a child cluster.
   *
   * Draw the masses \f$M'\f$ of the the cluster child produced 
   * by the fission of an heavy cluster (of mass M). m1, m2 are the masses
   * of the constituents of the cluster; m is the mass of the quark extract 
   * from the vacuum (together with its antiparticle). The algorithm produces
   * the mass of the cluster formed with consituent m1.
   * Two mass distributions can be used for the child cluster mass:
   * -# power-like mass distribution ("normal" mass) with power exp 
   *    \f[ M' = {\rm rnd}((M-m_1-m_2-m)^P, m^p)^{1/P} + m_1 \f]
   *    where \f$ P \f$ is a parameter of the model and \f$ \rm{rnd} \f$ is
   *    the function:
   *    \f[ \rm{rnd}(a,b) = (1-r)a + r b \f]
   *    and here \f$ r \f$ is a random number [0,1].
   * -# fast-decreasing exponential mass distribution ("soft" mass) with 
   *    rmin. rmin is given by 
   *    \f[ r_{\rm min} = \exp(-b (M - m_1 - m_2 - 2 m))  \f]
   *    where \f$ b \f$ is a parameter of the model. The generated mass is
   *    given by
   *    \f[ M' = m_1 + m - \frac{\log\left(
   *             {\rm rnd}(r_{\rm min}, 1-r_{\rm min})\right)}{b} \f].
   *
   * The choice of which mass distribution should be used for each of the two
   * cluster children is dictated by the parameter soft.
   */
  Energy drawChildMass(const Energy M, const Energy m1, const Energy m2, 
		       const Energy m, const double exp, const bool soft) const;

  /**
   * Determines the kinematics of a heavy cluster decay C->C1 + C2
   */
  void calculateKinematics(const Lorentz5Momentum &pClu, 
			   const Lorentz5Momentum &p0Q1, 
			   const bool toHadron1, const bool toHadron2,
			   Lorentz5Momentum &pClu1, Lorentz5Momentum &pClu2, 
			   Lorentz5Momentum &pQ1, Lorentz5Momentum &pQb, 
			   Lorentz5Momentum &pQ2, Lorentz5Momentum &pQ2b) const;

  /**
   * Determine the positions of the two children clusters.
   *
   * This routine generates the momentum of the decay products. It also
   * generates the momentum in the lab frame of the partons drawn out of
   * the vacuum.
   */
  void calculatePositions(const Lorentz5Momentum &pClu, 
		          const LorentzPoint & positionClu,
			  const Lorentz5Momentum & pClu1, 
			  const Lorentz5Momentum & pClu2, 
			  LorentzPoint & positionClu1, 
			  LorentzPoint & positionClu2 ) const;

protected:

  /** @name Access members for child classes. */
  //@{
  /**
   *  Access to the hadron selector
   */
  HadronSelectorPtr hadronsSelector() const {return _hadronsSelector;}

  /**
   *  Access to soft-cluster parameter
   */
  Energy btClM() const {return _btClM;}

  /**
   *  Cluster splitting paramater for light quarks
   */
  double pSplitLight() const {return _pSplitLight;}  

  /**
   *  Cluster splitting paramater for bottom quarks
   */
  double pSplitBottom() const {return _pSplitBottom;}

  /**
   *  Cluster splitting paramater for charm quarks
   */
  double pSplitCharm() const {return _pSplitCharm;}

  /**
   *  Cluster splitting paramater for exotic particles
   */
  double pSplitExotic() const {return _pSplitExotic;}
  //@}
  
private:

  /**
   * Check if a cluster is heavy enough to split again
   */
  bool isHeavy(tcClusterPtr );

  /**
   * A pointer to a Herwig::HadronSelector object for generating hadrons.
   */
  HadronSelectorPtr _hadronsSelector;

  /**
   * @name The Cluster max mass,dependant on which quarks are involved, used to determine when 
   * fission will occur.
   */
  //@{
  Energy _clMaxLight;
  Energy _clMaxBottom;
  Energy _clMaxCharm;
  Energy _clMaxExotic;
  //@}
  /**
   * @name The power used to determine when cluster fission will occur.
   */
  //@{
  double _clPowLight;
  double _clPowBottom;
  double _clPowCharm;
  double _clPowExotic;
  //@}
  /**
   * @name The power, dependant on whic quarks are involved, used in the cluster mass generation.
   */
  //@{
  double _pSplitLight;
  double _pSplitBottom;
  double _pSplitCharm;
  double _pSplitExotic;
  //@}
   /**
   * Parameter used (2/b) for the beam cluster mass generation.
   * Currently hard coded value.
   */
  Energy _btClM; 

  /**
   * Flag used to determine what distributions to use for the cluster masses.
   */
  int _iopRem;

  /**
   * The string constant
   */
  Tension _kappa;

};

}

#endif /* HERWIG_ClusterFissioner_H */
