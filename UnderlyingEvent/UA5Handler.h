// -*- C++ -*-
//
// UA5Handler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_UA5_H_
#define HERWIG_UA5_H_

#include <ThePEG/Handlers/HadronizationHandler.h>
#include "Herwig/Hadronization/CluHadConfig.h"
#include <ThePEG/Vectors/LorentzRotation.h>

namespace Herwig {

using namespace ThePEG;

/** \ingroup UnderlyingEvent
 *
 *  This is the class definition for the UA5Handler. This 
 *  class is designed to generate an underlying event     
 *  based on the UA5 model. This is intended as a basic   
 *  underlying event model which will be superceded by a  
 *  new model in Herwig.                                
 *                                                        
 *  This class interfaces with 
 *  the cluster hadronization. To that end there is an    
 *  interface set up with the ClusterFissioner class and  
 *  with the ClusterDecayer class.
 * 
 *  The Hadronization is responsible
 *  for the formation of the beam clusters. In this step the colour connection
 *  between the spectators and the initial-state parton showers is cut by the
 *  forced emission of a soft quark-antiquark pair. The underlying soft event in a
 *  hard hadron-hadron collision is then assumed to be a soft collision
 *  between these two beam clusters.
 *
 *  The model used for the underlying event is based on the minimum-bias
 *  \f$p\bar{p}\f$ event generator of the UA5 Collaboration,
 *  UA5 Collaboration, G.J. Alner et al., Nucl. Phys. B291 (1987) 445,
 *  modified to make use of our cluster fragmentation algorithm.
 *
 *  The parameter ProbSoft enables one to produce an underlying event in
 *  only a fraction ProbSoft of events (default=1.0).
 *
 *  The UA5 model starts from a parametrization of the \f$p\bar{p}\f$
 *  inelastic charged multiplicity distribution as a negative binomial distribution,
 * \f[P(n) = \frac{\Gamma(n+k)}{n!\,\Gamma(k)}
 *           \frac{(\bar n/k)^n}{(1+\bar n/k)^{n+k}}\;.\f]
 *  The parameter \f$k\f$ is given by
 * \f[1/k =k_1\ln s+k_2,\f]
 * and \f$\bar n\f$ is given by
 * \f[\bar n =n_1s^{n_2}+n_3\f]
 * As an option, for underlying events the value of \f$\sqrt{s}\f$ used to choose
 * the multiplicity \f$n\f$ may be increased by a factor EnhanceCM to allow
 * for an enhanced underlying activity in hard events.
 *
 * Once a charged multiplicity has been selected from the above distribution,
 * `softclusters' are generated with flavours \f$(f_1,f_2) = (q_{n-1},\bar q_n)\f$
 * by drawing \f$q_n\bar q_n = u\bar u$ or $d\bar d\f$ randomly from the 
 * vacuum. Soft cluster masses are assigned as
 * \f[M = m_{q1}+m_{q2}+m_1-\log(r_1 r_2)/m_2 \f]
 * where \f$r_{1,2}\f$ are random numbers, which gives a (shifted) exponential
 * distribution of \f$M^2\f$.  The parameters \f$m_1\f$ and \f$m_2\f$ control
 * the distribution and \f$m_{q1,2}\f$ are the masses of the quarks in the cluster.
 *
 * As each soft cluster is generated, it is decayed to stable hadrons using
 * the cluster hadronization model (without cluster fission) and the accumulated
 * charged multiplicity is computed.
 * Once the preselected charged multiplicity is exactly reached,
 * cluster generation is stopped. If it is exceeded, all clusters are rejected
 * and new ones are generated until the exact value is reached. In this way the
 * multiplicity distribution of stable charged hadrons
 * is generated exactly as prescribed.
 *
 * At this stage (to save time) the kinematic distribution of the soft clusters
 * has not yet been generated.  The decay products of each cluster are stored
 * in its rest frame.  Now the transverse momenta of the clusters are
 * generated with the distribution
 * \f[P(p_t)\propto p_t\exp\left(-p_{1,2,3}\sqrt{p_t^2+M^2}\right)\f]
 * where the slope parameter \f$p_{1,2,3}\f$ depends as indicated on the
 * flavour of the quark or diquark pair created in the
 * primary cluster decay, \f$p_1\f$ for light quarks, \f$p_2\f$ for the strange and
 * charm quarks and \f$p_3\f$ for diquarks.
 * Next the clusters are given a flat
 * rapidity distribution with Gaussian shoulders. The `reduced
 * rapidities' \f$\xi_i\f$ are generated first by drawing
 * from a distribution
 * \f[P(\xi) = N\;\;\;\mbox{for}\;|\xi|<0.6\f]
 * \f[P(\xi) = N\,e^{-(\xi-0.6)^2/2} \;\;\;\mbox{for}\;\xi>0.6\f]
 * \f[P(\xi) = N\,e^{-(\xi+0.6)^2/2} \;\;\;\mbox{for}\;\xi<-0.6\f]
 * where \f$N=1/(1.2+\sqrt{2\pi})\f$ is the normalization.  Next
 * a scaling factor \f$Y\f$ is computed such that the scaled cluster
 * rapidities \f$y_i=Y\xi_i\f$, their masses and transverse momenta
 * satisfy momentum conservation when compared to the total
 * energy of the underlying event. Thus the soft cluster rapidity
 * distribution retains its overall shape but becomes higher and
 * wider as the energy of the underlying event increases.
 *
 * Finally the decay products of each cluster are boosted from
 * its rest frame into the lab frame and added to the event record.
 *
 * @see \ref UA5HandlerInterfaces "The interfaces"
 * defined for UA5Handler..
 */

class UA5Handler : public HadronizationHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
   UA5Handler();
   //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

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

public:

  /**
   * This is the routine that starts the algorithm.
   * @param eh the EventHandler in charge of the generation.
   * @param tagged the vector of particles to consider. If empty, all
   * final state particles in the current Step is considered.
   * @param hint a possible Hint which is ignored in this implementation
   */
  virtual void handle(EventHandler &eh, const tPVector &tagged,
		      const Hint &hint) 
   ;

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
   *  Members to decay the clusters and hadrons produced in their decay,
   *  and insert the output in the event record.
   */
  //@{
  /**
   * Perform the decay of an unstable hadron.
   * @param parent the decaying particle
   * @param totalcharge The totalcharge of the decay proiducts
   * @param numbercharge The number of stabel charged decay products
   */
  void performDecay(PPtr parent,int & totalcharge,int & numbercharge) const;

  /**
   * Decay a cluster to two hadrons is sufficiently massive and to one if
   * not.
   * @param cluster The cluster to decay
   * @param single Whether or not ot allow decays to 
   */
  void decayCluster(ClusterPtr cluster, bool single) const;

  /**
   * Recursively add particle and decay products to the step
   * @param particle The particle
   * @param step The step
   * @param all Insert this particle as well as children
   */
  void insertParticle(PPtr particle,StepPtr step,bool all) const;
  //@}

  /**
   *  Members to generate the multiplicity according to a negative binomial
   *  distribution.
   */
  //@{
  /**
   * Calculate the negative binomial probability given the
   * mean \f$\bar n\f$, the multiplicity and \f$1/k\f$.
   * @param N  The multplicity for which to calculate the probability
   * @param mean The mean multiplicity \f$\bar n\f$
   * @param ek \f$1/k\f$
   * @return a value distributed according the negative binomial distribution
   */
  inline double negativeBinomial(int N, double mean, double ek) const;
  
  /**
   * The value of the mean multiplicity for a given energy E.
   * This is \f$n_1E^{2n_2}+n_3\f$ wher \f$n_1\f$, \f$n_2\f$ and \f$n_3\f$
   * are parameters.
   * @param E the energy to calculate the mean multiplicity for
   * @return the mean multiplicity
   */
  inline double meanMultiplicity(Energy E) const;
  
  /**
   * Generates a multiplicity for the energy E according to the negative
   * binomial distribution.
    * @param E The energy to generate for
    * @return the randomly generated multiplicity for the energy given
    */
  unsigned int multiplicity(Energy E) const;
  //@}
  
  /**
   * Members to generate the momenta of the clusters
   */
  //@{
  /**
   * This generates the momentum of the produced particles according to
   * the cylindrical phase space algorithm given
   * in Computer Physics Communications 9 (1975) 297-304 by S. Jadach.
   * @param clu1 The first incoming cluster
   * @param clu2 The second incoming cluster
   * @param clusters The list of clusters produced
   * @param CME The center of mass energy
   * @param cm The center of mass momentum (of the underlying event)
   */
  void generateMomentum(tClusterPtr clu1,tClusterPtr clu2,
			const ClusterVector &clusters, Energy CME,
			const Lorentz5Momentum & cm) const;
  
  /**
   * The implementation of the cylindrical phase space.
   * @param clusters The list of clusters to generate the momentum for
   * @param CME The center of mass energy
   */
  void generateCylindricalPS(const ClusterVector &clusters, Energy CME) const;
  //@}
  
  /**
   * This returns the rotation matrix needed to rotate p into the z axis
   */
  LorentzRotation rotate(const LorentzMomentum &p) const;
  
  /**
   *  Various methods to generate random distributions
   */
  //@{
  /**
   * Gaussian distribution
   * @param mean the mean of the distribution
   * @param stdev the standard deviation of the distribution
   * @return Arandom value from the gaussian distribution
   */
  template <typename T>
  inline T gaussDistribution(T mean, T stdev) const;
  
  /**
   * This returns a random number with a flat distribution
   * [-A,A] plus gaussian tail with stdev B 
   * TODO: Should move this to Utilities
   * @param A The width of the flat part
   * @param B The standard deviation of the gaussian tail
   * @return the randomly generated value
   */
  inline double randUng(double A, double B) const;
  
  /**
   * Generates a random azimuthal angle and puts x onto px and py 
   * TODO: Should move this to Utilities
   * @param pt The magnitude of the transverse momentum
   * @param px The x component after random rotation
   * @param py The y component after random rotation
   */
  template <typename T>
  inline void randAzm(T pt, T &px, T &py) const;
  
  /**
   * This returns random number from \f$dN/dp_T^2=exp(-p_{1,2,3}m_T\f$ distribution,
   * where \f$m_T=\sqrt{p_T^2+M^2}\f$. It uses Newton's method to solve \f$F-R=0\f$
   * @param AM0 The mass \f$M\f$.
   * @param B The slope
   * @return the value distributed from \f$dN/dp_T^2=exp(-p_{1,2,3}m_T\f$ with mean av
   */
  inline Energy randExt(Energy AM0,InvEnergy B) const;
  //@}
  
private:
  
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<UA5Handler> initUA5Handler;
  
  /**
   * This is never defined and since it can never be called it isn't 
   * needed. The prototype is defined so the compiler doesn't use the 
   * default = operator.
   */
  UA5Handler& operator=(const UA5Handler &);
  
private:

  /**
   *  Reference to the ClusterFissioner object
   */
  ClusterFissionerPtr clusterFissioner;
  
  /**
   *  Reference to the cluster decayer object.
   */
  ClusterDecayerPtr   clusterDecayer;
  
  /**
   *  Parameters for the mean multiplicity \f$\bar n =n_1s^{n_2}+n_3\f$
   */
  //@{
  /**
   *  The parameter \f$n_1\f$ 
   */
  double _n1;
  
  /**
   *  The parameter \f$n_2\f$ 
   */
  double _n2;
  
  /**
   *  The parameter \f$n_3\f$ 
   */
  double _n3;
  //@}
  
  /**
   *  Parameters for \f$k\f$ in the negative binomial 
   * distribution given by \f$1/k =k_1\ln s+k_2\f$
   */
  //@{
  /**
   *  The parameter \f$k_1\f$ 
   */
  double _k1;
  
  /**
   *  The parameter \f$k_2\f$ 
   */
  double _k2;
  //@}
  
  /**
   *  Parameters for the cluster mass distribution, 
   * \f$M = m_{q1}+m_{q2}+m_1-\log(r_1 r_2)/m_2\f$.
   */
  //@{
  /**
   *  The parameter \f$m_1\f$ 
   */
  Energy _m1;
  
  /**
   *  The parameter \f$m_2\f$ 
   */
  InvEnergy _m2;
  //@}
  
  /**
   *  Parameters for the transverpse momentum of the soft distribution,
   *  \f$P(p_T) \propto p_T*exp(-p_i\sqrt{p_T^2+M^2}\f$ 
   */
  //@{
  /**
   *  The parameter \f$p_1\f$ for light quarks
   */
  InvEnergy _p1;
  
  /**
   *  The parameter \f$p_2\f$ for strange and charm quarks
   */
  InvEnergy _p2;
  
  /**
   *  The parameter \f$p_3\f$ for diquarks
   */
  InvEnergy _p3;
  //@}
  
  /**
   * This is the probability of having a soft underlying event.
   */
  double _probSoft;
  
  /**
   * This is a parameter used to enhance the CM energy used to 
   * generate the multiplicity distribution.
   */
  double _enhanceCM; 
  
  /**
   *  The maximum number of attempts to generate the distribution
   */
  unsigned int _maxtries;
  

  /**
   * Whether to warn about using UA5 and MPI simultaneously.
   */
  bool _needWarning;
};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of UA5Handler. */
template<>
struct BaseClassTrait<Herwig::UA5Handler,1> { 
  /** Typedef of the first base class of UA5Handler. */
  typedef HadronizationHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the UA5Handler class and the shared object where it is defined. */
template<>
struct ClassTraits<Herwig::UA5Handler> :
  public ClassTraitsBase<Herwig::UA5Handler> {
  /** Return a platform-independent class name */
    static string className() { return "Herwig::UA5Handler"; }
  /** Return the name of the shared library be loaded to get
   *  access to the WeakPartonicDecayer class and every other class it uses
   *  (except the base class). */
    static string library() { return "HwUA5.so"; }
};

/** @endcond */

}

#include "UA5Handler.icc"

#endif
