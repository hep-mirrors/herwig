// -*- C++ -*-
#ifndef HERWIG_ClusterFissioner_H
#define HERWIG_ClusterFissioner_H
/*! \class Herwig::ClusterFissioner ClusterFissioner.h "Herwig++\Hadronization\ClusterFissioner.h"
 * \brief This class handles clusters which are too heavy.
 * \author Philip Stephens
 * \author Alberto Ribon
 * \author Stefan Gieseke
 * \ingroup Hadronization
 *
 * This class does the job of chopping up either heavy clusters or beam 
 * clusters in two lighter ones. The procedure is repeated recursively until 
 * all of the cluster children have masses below some threshold values.
 *
 * For the beam remnant clusters, at the moment what is done is the following.
 * In the case that the soft underlying event is switched on, the 
 * beam remnant clusters are tagged as not available,
 * therefore they will not be treated at all during the hadronization. 
 * In the case instead that the soft underlying event is switched off,
 * then the beam remnant clusters are treated exactly as "normal" clusters,
 * with the only exception of the mass spectrum used to generate the
 * cluster children masses. For non-beam clusters, the masses of the cluster
 * children are draw from a power-like mass distribution; for beam clusters,
 * according to the value of the flag _IOpRem, either both 
 * children masses are draw from a fast-decreasing exponential mass 
 * distribution (case _IOpRem == 0, or, indendently by 
 * _IOpRem, in the special case that the beam cluster contains two 
 * beam remnants), or one mass from the exponential distribution (corresponding
 *  of the cluster child with the beam remnant) and the other with the usual 
 * power-like distribution (case _IOpRem == 1, which is the 
 * default one, as in Herwig 6.3). 
 *
 * The reason behind the use of a fast-decreasing exponential distribution 
 * is that to avoid a large transverse energy from the many sequential
 * fissions that would otherwise occur due to the typical large cluster 
 * mass of beam clusters. Using instead an exponential distribution 
 * the masses of the two cluster children will be very small (order of 
 * GeV).
 *
 * The rationale behind the implementation of the splitting of clusters
 * has been to preserve *all* of the information about such splitting 
 * process. More explicitly a ThePEG::Step class is passed in and the
 * new clusters are added to the step as the decay products of the
 * heavy cluster. This approach has the twofold 
 * advantage to provide all of the information that could be needed 
 * (expecially in future developments), without any information loss, 
 * and furthermore it allows a better debugging. 
 *
 * See also:
 * GlobalParameters.h, HadronSelector.h.
 */ 

#include <ThePEG/Handlers/HandlerBase.h>
#include "CluHadConfig.h"
#include "HadronSelector.h"
#include "Herwig++/Utilities/GlobalParameters.h"


namespace Herwig {


using namespace ThePEG;

  //class Cluster;          // forward declaration


class ClusterFissioner: public ThePEG::HandlerBase {

public:

  inline ClusterFissioner();
  inline ClusterFissioner(const ClusterFissioner &);
  virtual ~ClusterFissioner();
  // Standard ctors and dtor.

  void fission(const StepPtr &);
  /*!< Splits the clusters which are too heavy.
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
    
public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  static void Init();
  //!< Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  //!< Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  //!< Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<ClusterFissioner> initClusterFissioner;
  //!< Describe a concrete class with persistent data.

  ClusterFissioner & operator=(const ClusterFissioner &);
  //!< Private and non-existent assignment operator.

  void cut(tClusterPtr, const StepPtr&, ClusterVector&);
  /*!< This method directs the splitting of the heavy clusters
   *
   * This method does the splitting of the clusters and all of its cluster 
   * children, if heavy. All of these new children clusters are added to the
   * collection of clusters. The method works as follows.
   * Initially the vector contains just the input pointer to the
   * cluster to be split. Then it will be filled recursively by all
   * of the cluster's children that are heavy enough to require
   * to be split. In each loop, the last element of the vector is 
   * considered (only once because it is then removed from the vector).
   * This approach is conceptually recursive, but avoids the overhead of
   * a concrete recursive function. Furthermore it requires minimal changes
   * in the case that the fission of an heavy cluster could produce more
   * than two cluster children as assumed now. 
   *
   * For normal, non-beam clusters, a power-like mass distribution
   * is used, whereas for beam clusters a fast-decreasing exponential mass 
   * distribution is used instead. This avoids many iterative splitting which 
   * could produce an unphysical large transverse energy from a supposed 
   * soft beam remnant process.
   */
public:
  typedef pair<PPtr,PPtr> PPair;
  //!< Definition for easy passing of two particles.
  typedef pair<PPair,PPair> cutType;
  //!< Definition for use in the cut function.
  cutType cut(tClusterPtr &);
  /*!< Splits the input cluster.
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

private:
  PPair produceHadron(const long id1, const long id2, Lorentz5Momentum &a,
		      LorentzPoint &b) const;
  /*!< Produces a hadron and returns the flavour drawn from the vacuum.
   *
   * This routine produces a new hadron from the flavours id1 and id2. It
   * also sets the momentum and vertex to the values given.
   */
  PPair produceCluster(tPPtr &p1, const long id, Lorentz5Momentum &a, 
		       LorentzPoint &b, Lorentz5Momentum &c, 
		       Lorentz5Momentum &d, const bool rem) const;
  /*!< Produces a cluster from the flavours passed in.
   *
   * This routine produces a new cluster with the flavours given by p1 and id.
   * The new 5 momentum is a and the parent momentum are c and d. C is for the
   * p1 and d is for the new particle id. rem specifies whether the existing
   * particle is a beam remnant or not.
   */

  long drawNewFlavour() const;
  /*!< Returns a flavour from the choice u,d,s.
   *
   * Return the id ( >0 ) of the quark-antiquark pair from the vacuum
   * needed for fission of a heavy cluster. Equal probabilities
   * are assumed for quarks  u , d , s . 
   */

  void drawChildMass(const Energy M, const Energy m1, const Energy m2, 
		     const Energy m, Energy & Mc, const double exp,
                           const double b, const bool rem) const; 
  /*!< Produces the mass of a child cluster.
   *
   * Draw the masses M' of the the cluster child produced 
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
   * cluster children is dictated by the bool rem. If _IOpRem is 0, the
   * soft distribution is always used.
   *
   * Finally, sometimes, when the phase space available is tiny, many attempts 
   * fail to produce a pair of masses kinematically acceptable; in these cases 
   * it gives up returning false, otherwise it returns true when the splitting 
   * succeeds.
   */

  void calculateKinematics(const Lorentz5Momentum &pClu, 
		           const Lorentz5Momentum &p0Q1, 
			   const bool toHadron1, const bool toHadron2,
			   Lorentz5Momentum &pClu1, Lorentz5Momentum &pClu2, 
			   Lorentz5Momentum &pQ1, Lorentz5Momentum &pQb, 
			   Lorentz5Momentum &pQ2, Lorentz5Momentum &pQ2b) const;
  //!< Determines the kinematics of a heavy cluster decay C->C1 + C2

  
  void calculatePositions(const Lorentz5Momentum &pClu, 
		          const LorentzPoint & positionClu,
			  const Lorentz5Momentum & pClu1, 
			  const Lorentz5Momentum & pClu2, 
			  LorentzPoint & positionClu1, 
			  LorentzPoint & positionClu2 ) const;
  /*!< Determine the positions of the two children clusters.
   *
   * This routine generates the momentum of the decay products. It also
   * generates the momentum in the lab frame of the partons drawn out of
   * the vacuum.
   */

  HadronSelectorPtr _hadronsSelector;
  //!< A pointer to a Herwig::HadronSelector object for generating hadrons.
  GlobParamPtr      _globalParameters;  
  //!< A pointer to a Herwig::GlobalParameters object for global variables.

  Energy _ClMax;
  //!< The Cluster max mass used to determine when fission will occur.
  double _ClPow;
  //!< The power used to determine when cluster fission will occur.
  double _PSplt1;
  //!< The power used in the cluster mass generation. This is the non-b param.
  double _PSplt2;
  //!< The power used in the cluster mass generation. This is the b param.

  Energy _BtClM; 
  //!< Parameter used (2/b) for the beam cluster mass generation. Currently hard coded value.
  int _IOpRem;
  //!< Flag used to determine what distributions to use for the cluster masses.

};

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of ClusterFissioner.
template <>
struct BaseClassTrait<Herwig::ClusterFissioner,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ClusterFissioner>: public ClassTraitsBase<Herwig::ClusterFissioner> {
  static string className() { return "/Herwig++/ClusterFissioner"; }
  // Return the class name.
  static string library() { return "libHwHadronization.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#endif // DOXYGEN

#include "ClusterFissioner.icc"

#endif /* HERWIG_ClusterFissioner_H */
