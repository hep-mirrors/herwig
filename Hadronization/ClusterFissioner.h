// -*- C++ -*-
//
// ClusterFissioner.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ClusterFissioner_H
#define HERWIG_ClusterFissioner_H

#include <ThePEG/Interface/Interfaced.h>
#include "CluHadConfig.h"
#include "ClusterFissioner.fh"
#include "HadronSpectrum.h"

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
  /**
   * Default destructor.
   */
   virtual ~ClusterFissioner();

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
  virtual tPVector fission(ClusterVector & clusters, bool softUEisOn);

  /**
   * Return the hadron spectrum
   */
  virtual Ptr<HadronSpectrum>::tptr spectrum() const {
    return _hadronSpectrum;
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

private:

  /**
   * Private and non-existent assignment operator.
   */
  ClusterFissioner & operator=(const ClusterFissioner &) = delete;

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
  void drawNewFlavourQuarks(PPtr& newPtrPos,PPtr& newPtrNeg) const;

  /**
   * Returns the new quark-antiquark pair or diquark -
   * antidiquark pair needed for fission of a heavy cluster.
   */
  void drawNewFlavourDiquarks(PPtr& newPtrPos,PPtr& newPtrNeg,
		  		const ClusterPtr & clu) const;

  /**
   * Returns the new quark-antiquark pair
   * needed for fission of a heavy cluster. Equal probabilities
   * are assumed for producing  u, d, or s pairs.
   * Extra argument is used when performing strangeness enhancement
   */
  void drawNewFlavourEnhanced(PPtr& newPtrPos,PPtr& newPtrNeg, Energy2 mass2) const;


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

  /**
   *  Dimension used to calculate phase space weights
   */
  double dim() const {return _dim;}

  /**
   *  Access to soft-cluster parameter
   */
  Energy btClM() const {return _btClM;}

  /**
  *  Function that returns either the cluster mass or the lambda measure
  */
  Energy2 clustermass(const ClusterPtr & cluster) const;
  
  /**
   * Draw a new flavour for the given cluster; currently defaults to
   * the default model
   */
  virtual void drawNewFlavour(PPtr& newPtr1, PPtr& newPtr2, const ClusterPtr & cluster) const {
    if (_enhanceSProb == 0){
      if (_diquarkClusterFission>=0) drawNewFlavourDiquarks(newPtr1,newPtr2,cluster);
			else drawNewFlavourQuarks(newPtr1,newPtr2);
    }
    else {
      drawNewFlavourEnhanced(newPtr1,newPtr2,clustermass(cluster));
    }
  }

  /**
   * Calculate the masses and possibly kinematics of the cluster
   * fission at hand; if claculateKineamtics is perfomring non-trivial
   * steps kinematics claulcated here will be overriden. Currentl;y resorts to the default
	 * @return the potentially non-trivial distribution weight=f(M1,M2)
	 *         On Failure we return 0
   */
  virtual double drawNewMasses(const Energy Mc, const bool soft1, const bool soft2,
					    Lorentz5Momentum& pClu1, Lorentz5Momentum& pClu2,
					    tcPPtr ptrQ1,    const Lorentz5Momentum& pQ1, 
					    tcPPtr, const Lorentz5Momentum& pQone,
					    tcPPtr, const Lorentz5Momentum& pQtwo,
					    tcPPtr ptrQ2,    const Lorentz5Momentum& pQ2) const;

  /**
   * Calculate the final kinematics of a heavy cluster decay C->C1 +
   * C2, if not already performed by drawNewMasses
   */
	virtual void calculateKinematics(const Lorentz5Momentum & pClu,
					   const Lorentz5Momentum & p0Q1,
					   const bool toHadron1,
					   const bool toHadron2,
					   Lorentz5Momentum & pClu1,
					   Lorentz5Momentum & pClu2,
					   Lorentz5Momentum & pQ1,
					   Lorentz5Momentum & pQbar,
					   Lorentz5Momentum & pQ,
					   Lorentz5Momentum & pQ2bar) const;

protected:

  /** @name Access members for child classes. */
  //@{
  /**
   *  Access to the hadron selector
   */
  HadronSpectrumPtr hadronSpectrum() const {return _hadronSpectrum;}
  /**
   *  Access for fission Pwts
   */
  const map<long,double> fissionPwt() const { return _fissionPwt;}

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
	 * Flat PhaseSpace weight for ClusterFission
	 */
	double weightFlatPhaseSpace(const Energy Mc, const Energy Mc1, const Energy Mc2,
			const Energy m, const Energy m1, const Energy m2,
			tcPPtr pQ, tcPPtr pQ1, tcPPtr pQ2) const {
		switch (_phaseSpaceWeights)
		{
			case 1:
				return weightPhaseSpaceConstituentMasses(Mc, Mc1, Mc2, m, m1, m2, 0.0);
			case 2:
				return weightFlatPhaseSpaceHadronMasses(Mc, Mc1, Mc2, pQ, pQ1, pQ2);
			case 3:
				return weightFlatPhaseSpaceNoConstituentMasses(Mc, Mc1, Mc2);
			default:
				assert(false);
		}
	};
	/**
	 * PhaseSpace weight for ClusterFission using constituent masses
	 */
	double weightPhaseSpaceConstituentMasses(const Energy Mc, const Energy Mc1, const Energy Mc2,
			const Energy m, const Energy m1, const Energy m2, const double power=0.0) const;
	/**
	 * Flat PhaseSpace weight for ClusterFission using lightest hadron masses
	 */
	double weightFlatPhaseSpaceHadronMasses(const Energy Mc, const Energy Mc1, const Energy Mc2,
			tcPPtr pQ, tcPPtr pQ1, tcPPtr pQ2) const;
	double weightFlatPhaseSpaceNoConstituentMasses(const Energy Mc, const Energy Mc1, const Energy Mc2) const;


  /**
   * Calculate a veto for the phase space weight
   */
  bool phaseSpaceVeto(const Energy Mc, const Energy Mc1, const Energy Mc2,
			     const Energy m, const Energy m1, const Energy m2, tcPPtr pQ1=tcPPtr(), tcPPtr pQ2=tcPPtr(), tcPPtr pQ=tcPPtr(), const double power = 0.0) const;
  bool phaseSpaceVetoConstituentMasses(const Energy Mc, const Energy Mc1, const Energy Mc2,
			     const Energy m, const Energy m1, const Energy m2, const double power = 0.0) const;
  bool phaseSpaceVetoNoConstituentMasses(const Energy Mc, const Energy Mc1, const Energy Mc2) const;

  bool phaseSpaceVetoHadronPairs(const Energy Mc, const Energy Mc1, const Energy Mc2,
			      tcPPtr pQ1, tcPPtr pQ2, tcPPtr pQconst) const;


//@}

protected:

  /**
  * Smooth probability for dynamic threshold cuts:
  * @scale the current scale, e.g. the mass of the cluster,
  * @threshold the physical threshold,
   */
  bool ProbabilityFunction(double scale, double threshold);
	bool ProbabilityFunctionPower(double Mass, double threshold);

  /**
   * Check if a cluster is heavy enough to split again
   */
  bool isHeavy(tcClusterPtr );

  /**
   * Check if a cluster is heavy enough to be at least kinematically able to split
   */
  bool canSplitMinimally(tcClusterPtr, Energy);

  /**
   *  Check if can't make a hadron from the partons
   */
  inline bool cantMakeHadron(tcPPtr p1, tcPPtr p2) {
    return ! spectrum()->canBeHadron(p1->dataPtr(), p2->dataPtr());
  }
  /**
   * A pointer to a Herwig::HadronSpectrum object for generating hadrons.
   */
  HadronSpectrumPtr _hadronSpectrum;

  /**
   * @name The Cluster max mass,dependant on which quarks are involved, used to determine when
   * fission will occur.
   */
  //@{
  Energy _clMaxLight;
  Energy _clMaxDiquark;
  map<long,Energy> _clMaxHeavy;
  Energy _clMaxExotic;
  //@}
  /**
   * @name The power used to determine when cluster fission will occur.
   */
  //@{
  double _clPowLight;
  double _clPowDiquark;
  map<long,double> _clPowHeavy;
  double _clPowExotic;
  //@}
  /**
   * @name The power, dependant on whic quarks are involved, used in the cluster mass generation.
   */
  //@{
  double _pSplitLight;
  map<long,double> _pSplitHeavy;
  double _pSplitExotic;

  /**
   * Weights for alternative cluster fission
   */
  map<long,double> _fissionPwt;

  /**
   * Include phase space weights
   */
  int _phaseSpaceWeights;

  /**
   * Dimensionality of phase space weight
   */
  double _dim;

  /**
  * Flag used to determine between normal cluster fission and alternative cluster fission
  */
  int _fissionCluster;

  /**
  * Flag to choose static or dynamic kinematic thresholds in cluster splittings
  */
  int _kinematicThresholdChoice;

  /**
   * Pwt weight for drawing diquark
   */
  double _pwtDIquark;

  /**
   * allow clusters to fission to 1 (or 2) diquark clusters or not
   */
  int _diquarkClusterFission;

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

  /**
  *  Flag that switches between no strangeness enhancement, scaling enhancement,
  *  and exponential enhancement (in numerical order)
  */
  int _enhanceSProb;

  /**
  *  Parameter that governs the strangeness enhancement scaling
  */
  Energy _m0Fission;

  /**
  *  Flag that switches between mass measures used in strangeness enhancement:
  *  cluster mass, or the lambda measure -  ( m_{clu}^2 - (m_q + m_{qbar})^2 )
  */
  int _massMeasure;

  /**
  *  Constant variable which stops the scale from being to large, and not worth
  *  calculating
  */
  const double _maxScale = 20.;

 /**
  * Power factor in ClausterFissioner bell probablity function
  */
  double _probPowFactor;

  /**
  * Shifts from the center in ClausterFissioner bell probablity function
  */
  double _probShift;

  /**
  * Shifts from the kinetic threshold in ClausterFissioner
  */
  Energy2 _kinThresholdShift;

  /**
   * Flag for strict diquark selection according to kinematics
   */
  int _strictDiquarkKinematics;

  /**
   * Use Covariant boost in MatrixElementClusterFissioner
   */
  bool _covariantBoost;

	/**
	 * Power for MassPreSampler = PowerLaw
	 */
	double _powerLawPower;
	
	/*
	 * flag for allowing strange Diquarks to be produced during
	 * Cluster Fission
	 * */
	unsigned int _hadronizingStrangeDiquarks;

	private:
	/*
	 * DEBUG output */
	int _writeOut;
};

}

#endif /* HERWIG_ClusterFissioner_H */
