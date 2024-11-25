// -*- C++ -*-
//
// MatrixElementClusterFissioner.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatrixElementClusterFissioner_H
#define HERWIG_MatrixElementClusterFissioner_H

#include <ThePEG/Interface/Interfaced.h>
#include "CluHadConfig.h"
#include "ClusterFissioner.h"
#include "HadronSpectrum.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Hadronization
 *  \class MatrixElementClusterFissioner
 *  \brief This class handles clusters which are too heavy by a semi-perturbative cluster fission matrix element.
 *  \author Stefan Kiebacher
 *
 * @see \ref MatrixElementClusterFissionerInterfaces "The interfaces"
 * defined for MatrixElementClusterFissioner.
 */
class MatrixElementClusterFissioner: public ClusterFissioner {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
   MatrixElementClusterFissioner();
  /**
   * Default destructor.
   */
   virtual ~MatrixElementClusterFissioner();

  //@}
  virtual tPVector fission(ClusterVector & clusters, bool softUEisOn);

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
  MatrixElementClusterFissioner & operator=(const MatrixElementClusterFissioner &) = delete;

  virtual void cut(stack<ClusterPtr> &,
	   ClusterVector&, tPVector & finalhadrons, bool softUEisOn);

public:

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
   *  Split two-component cluster using New FissionApproach
   */
  virtual cutType cutTwoMatrixElement(ClusterPtr &, tPVector & finalhadrons, bool softUEisOn);

  //@}

protected:
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
   * Calculate the masses and possibly kinematics of the cluster
   * fission at hand; if claculateKineamtics is perfomring non-trivial
   * steps kinematics claulcated here will be overriden. Currentl;y resorts to the default
   */
  double drawNewMassesDefault(const Energy Mc, const bool soft1, const bool soft2,
					    Lorentz5Momentum& pClu1, Lorentz5Momentum& pClu2,
					    tcPPtr ptrQ1, const Lorentz5Momentum& pQ1, 
					    tcPPtr, const Lorentz5Momentum& pQone,
					    tcPPtr, const Lorentz5Momentum& pQtwo,
					    tcPPtr ptrQ2,  const Lorentz5Momentum& pQ2) const;

  /**
   * Sample the masses for flat phase space
	 * */
  double drawNewMassesUniform(const Energy Mc,
					    Lorentz5Momentum & pClu1, Lorentz5Momentum & pClu2,
					    const Lorentz5Momentum & pQ1, 
					    const Lorentz5Momentum & pQ,
					    const Lorentz5Momentum & pQ2) const;

  /**
   * Sample the masses for flat phase space
	 * */
  double drawNewMassesPhaseSpace(const Energy Mc,
					    Lorentz5Momentum & pClu1, Lorentz5Momentum & pClu2,
					    const Lorentz5Momentum & pQ1, 
					    const Lorentz5Momentum & pQ,
					    const Lorentz5Momentum & pQ2,
							tcPPtr ptrQ1, tcPPtr ptrQ2, tcPPtr ptrQ ) const;
  /**
   * Sample the masses for flat phase space modulated by a power law
	 * */
	double drawNewMassesPhaseSpacePowerLaw(const Energy Mc,
			Lorentz5Momentum& pClu1, Lorentz5Momentum& pClu2,
			const Lorentz5Momentum& pQ1, 
			const Lorentz5Momentum& pQ,
			const Lorentz5Momentum& pQ2,
			tcPPtr ptrQ1, tcPPtr ptrQ2, tcPPtr ptrQ) const;

  /**
   * Calculate the final kinematics of a heavy cluster decay C->C1 +
   * C2, if not already performed by drawNewMasses
	 * @return  returns false if failes
   */
  double drawKinematics(const Lorentz5Momentum &pClu,
				   const Lorentz5Momentum &p0Q1,
				   const Lorentz5Momentum &p0Q2,
				   const bool toHadron1, const bool toHadron2,
				   Lorentz5Momentum &pClu1, Lorentz5Momentum &pClu2,
				   Lorentz5Momentum &pQ1, Lorentz5Momentum &pQb,
				   Lorentz5Momentum &pQ,  Lorentz5Momentum &pQ2b) const;

	/**
	 * Calculation of the squared matrix element from the drawn 
	 * momenta
	 * @return value of the squared matrix element in units of
	 *         1/GeV4
	 * */
	double calculateSQME(
			const Lorentz5Momentum & p1,
			const Lorentz5Momentum & p2,
			const Lorentz5Momentum & q1,
			const Lorentz5Momentum & q,
			const Lorentz5Momentum & q2,
			const Lorentz5Momentum & qbar) const;

	/**
	 * Calculation of the overestimate for the squared matrix 
	 * element independent on M1 and M2
	 * @return value of the overestimate squared matrix element
	 *         in units of 1/GeV4
	 * */
	double calculateSQME_OverEstimate(
			const Energy& Mc,
			const Energy& m1,
			const Energy& m2,
			const Energy& mq) const;
protected:

	/** @name Access members for child classes. */
	//@{
  /**
   *  Access to the hadron selector
   */
  // HadronSpectrumPtr hadronSpectrum() const {return _hadronSpectrum;}
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

	std::pair<Axis,double> sampleDirectionCluster(const Lorentz5Momentum & pQ, const Lorentz5Momentum & pClu)  const;
	std::pair<Axis,double> sampleDirectionConstituents(const Lorentz5Momentum & pClu, const Energy Mcluster)  const;
  /**
   * Samples the direction of Cluster Fission products uniformly
   **/
  Axis sampleDirectionIsotropic() const;

  /**
   * Samples the direction of Cluster Fission products uniformly
   * but only accepts those flying in the direction of pRelCOM
   **/
  Axis sampleDirectionSemiUniform(const Lorentz5Momentum & pQ, const Lorentz5Momentum & pClu) const;

  /**
   * Samples the direction of Cluster Fission products according to Default
   * fully aligned D = 1+1 Fission
   * */
  Axis sampleDirectionAligned(const Lorentz5Momentum & pQ, const Lorentz5Momentum & pClu) const;

  /**
   * Samples the direction of Cluster Fission products according to Default
   * Tchannel like distribution
   * */
  std::pair<Axis,double>  sampleDirectionTchannel(const Axis & dirQ, const double Ainv) const;
  
  /**
   * Samples the direction of direction from an exponential distribution
   * */
	std::pair<Axis,double> sampleDirectionExponential(const Axis & dirQ, const double lambda) const;

	/**
   * Samples the direction of Cluster Fission products according to smeared direction along the cluster
	 * */
	Axis sampleDirectionAlignedSmeared(const Lorentz5Momentum & pQ, const Lorentz5Momentum & pClu) const;

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
	
	void countPaccGreater1(){
		_counterPaccGreater1++;
	};
	void countMaxLoopViolations(){
		_counterMaxLoopViolations++;
	};
	void countFissionMatrixElement(){
		_counterFissionMatrixElement++;
	};

private:
  /**
   * Flag for allowing Hadron Final states in Cluster Fission
   */
  int _allowHadronFinalStates;

  /**
   * Choice of Mass sampling for MatrixElementClusterFissioner
	 * i.e. rejection sampling starting from MassPresampling
	 * Chooses the to be sampled Mass distribution
	 * Note: This ideal distribution would be ideally
	 *       exactly the integral over the angles as a
	 *       function of M1,M2
   */
  int _massSampler;

  /**
   * Choice of Phase Space sampling for MatrixElementClusterFissioner
	 * for child cluster directions and constituent directions
   */
  int _phaseSpaceSamplerCluster;
  int _phaseSpaceSamplerConstituent;

  /**
   * Choice of Matrix Element for MatrixElementClusterFissioner
	 * Note: The choice of the matrix element requires
	 *       to provide one overestimate as a function
	 *       of M1,M2 and another overestimate of the
	 *       previous overestimate independent of 
	 *       M1 and M2 i.e. only dependent on M,m1,m2,m
   */
  int _matrixElement;

  /**
   * Choice of MatrixElementClusterFissioner Approach
   */
	int _fissionApproach;

	/**
	 * Power for MassPreSampler = PowerLaw
	 */
	double _powerLawPower;

  /**
   * Choice of MatrixElementClusterFissioner Approach
	 * Technical Parameter for how many tries are allowed to sample the
	 * Cluster Fission matrix element before reverting to fissioning 
	 * using the default Fission Aproach
   */
	int _maxLoopFissionMatrixElement;

	/* 
	 * Safety factor for a better overestimate of the matrix Element
	 *
	 * */
	double _safetyFactorMatrixElement;

	/* 
	 * IR cutOff for IR divergent diagramms
	 * */
	double _epsilonResolution;
	
	/* 
	 * Failing mode for MaxLoopMatrixElement
	 * */
	int _failModeMaxLoopMatrixElement;
	// For MatrixElement only:

	// Counter for how many times we accept events with Pacc>1
	static unsigned int _counterPaccGreater1;
	// Counter for how many times we escape the sampling by
	// tries > _maxLoopFissionMatrixElement
	// NOTE: This maxing out either rejects event or uses old
	//       ClusterFission 
	static unsigned int _counterMaxLoopViolations;
	static unsigned int _counterFissionMatrixElement;
	private:
	int _writeOut;

};

}

#endif /* HERWIG_MatrixElementClusterFissioner_H */
