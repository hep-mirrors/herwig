// -*- C++ -*-
#ifndef HERWIG_ClusterFissioner_H
#define HERWIG_ClusterFissioner_H
//
// This is the declaration of the <!id>ClusterFissioner<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class does the job of chopping up either heavy clusters or beam 
// clusters in two lighter ones. The procedure is repeated recursively until 
// all of the cluster children have masses below some threshold values. <BR>
//
// For the beam remnant clusters, at the moment what is done is the following.
// In the case that the soft underlying event is switched on, the 
// beam remnant clusters are tagged as not available,
// therefore they will not be treated at all during the hadronization. 
// In the case instead that the soft underlying event is switched off,
// then the beam remnant clusters are treated exactly as "normal" clusters,
// with the only exception of the mass spectrum used to generate the
// cluster children masses. For non-beam clusters, the masses of the cluster
// children are drawn from a power-like mass distribution; for beam clusters,
// according to the value of the flag <!id>_IOpRem<!!id>, either both 
// children masses are drawn from a fast-decreasing exponential mass 
// distribution (case <!id>_IOpRem == 0<!!id>, or, indendently by 
// <!id>_IOpRem<!!id>, in the special case that the beam cluster contains two 
// beam remnants), or one mass from the exponential distribution (corresponding
//  of the cluster child with the beam remnant) and the other with the usual 
// power-like distribution (case <!id>_IOpRem == 1<!!id>, which is the 
// default one, as in Herwig 6.3).  <BR>
// The reason behind the use of a fast-decreasing exponential distribution 
// is that to avoid a large transverse energy from the many sequential
// fissions that would otherwise occur due to the typical large cluster 
// mass of beam clusters. Using instead an exponential distribution 
// the masses of the two cluster children will be very small (order of 
// <I>GeV</I>). <BR>
//
// The rationale behind the implementation of the splitting of clusters
// has been to preserve *all* of the information about such splitting 
// process. More explicitly, at the end of the full splitting, the 
// container of cluster (pointers) <!id>_collecCluPtr<!!id>  (which is passed 
// to this class <!class>ClusterFissioner<!!class> by the main one 
// <!class>ClusterHadronizationHandler<!!class>) 
// does not have only the final cluster products (the onces not heavy 
// and therefore that not need split, and that are ready to be decayed 
// in hadrons) but all of the clusters, including the initial ones, 
// formed during the cluster finding stage, and all of the intermediate 
// heavy clusters that have been split. This approach has the twofold 
// advantage to provide all of the information that could be needed 
// (expecially in future developments), without any information loss, 
// and furthermore it allows a better debugging. There is, however, a 
// small price to pay: in the following stage, that is the decay of 
// clusters into hadrons, we have to iterate over all of the clusters,
// and not directly to the final ones, although we strictly need to work 
// only on the latter. We think this is a minimal overhead. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:GlobalParameters.html">GlobalParameters.h</a>, <BR>
// <a href="http:HadronsSelector.html">HadronsSelector.h</a>.
// 

#include "Pythia7/Handlers/HandlerBase.h"
#include "CluHadConfig.h"
#include "HadronsSelector.h"
#include "Herwig++/Config/GlobalParameters.h"


namespace Herwig {


using namespace Pythia7;

  //class Cluster;          // forward declaration


class ClusterFissioner: public Pythia7::HandlerBase {

public:

  inline ClusterFissioner();
  inline ClusterFissioner(const ClusterFissioner &);
  virtual ~ClusterFissioner();
  // Standard ctors and dtor.

  void fission(const StepPtr &);
  // Split either heavy clusters or beam clusters recursively until all 
  // children have mass below some threshold. 
  // For beam clusters, they are split only if the soft underlying event
  // is switched off, otherwise these clusters will be tagged as unavailable
  // and they will not be treated by the hadronization altogether. 
  // In the case beam clusters will be split, the procedure is exactly
  // the same as for normal non-beam clusters, with the only exception
  // of the mass spectrum from which to draw the masses of the two 
  // cluster children (see method <!id>drawnChildrenMasses<!!id> for details).
    
public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<ClusterFissioner> initClusterFissioner;
  // Describe a concrete class with persistent data.

  ClusterFissioner & operator=(const ClusterFissioner &);
  // Private and non-existent assignment operator.

  void cut(tClusterPtr cluPtr, const StepPtr&, ClusterVector &clusters);
  // Split the input cluster (which can be either an heavy non-beam
  // cluster or a beam cluster) recursively until all children have
  // mass below some threshold. The "output" consists of all children and
  // grand-children clusters that are coming from the fission of the input 
  // cluster: all of these clusters (indeed the pointers to them) are added 
  // in the container  collecCluPtr. In the (rare) cases in which 
  // one of the two fission products is light (that is below the 
  // mass of the two lightest hadrons with proper flavour numbers)
  // then a hadron (the lightest one) is created instead of a cluster
  // and the parent heavy cluster has a single cluster child and a
  // single hadron child, rather than two cluster children as usual.
  // Notice that this method treats also beam clusters, besides
  // heavy non-beam clusters: the difference between these two types
  // of clusters is only in the mass spectrum used to generate the
  // cluster children masses (see method drawnChildrenMasses).

  long drawnNewFlavour() const;
  // Return the id ( >0 ) of the quark-antiquark pair from the vacuum
  // needed for fission of an heavy cluster. Equal probabilities
  // are assumed for quarks  u , d , s . 

  bool drawnChildrenMasses(const Energy Mclu, const Energy m1, const Energy m2, 
			   const Energy m3, Energy & Mclu1, Energy & Mclu2,
			   const double exponent1, const double exponent2,
                           const Energy average, const int iRemnant) const; 
  // Drawn the masses (Mclu1, Mclu2) of the the two clusters children produced 
  // by the fission of an heavy cluster (of mass Mclu). m1, m2 are the masses
  // of the constituents of the cluster; m3 is the mass of the quark extract from
  // the vacuum (together with its antiparticle). 
  // Two mass distributions can be used for the children cluster masses:
  // 1) power-like mass distribution ("normal" mass) with power exponent1 
  //    and exponent2 for the two cluster children respectively; 
  // 2) fast-decreasing exponential mass distribution ("soft" mass) with 
  //    average given by  average  method parameter.
  // The choice of which mass distribution should be used for each of the two
  // cluster children is dictated by the integer  iRemnant, as follows:
  // i)   iRemnant == 0  then both Mclu1 and Mclu2 are "normal" masses;
  // ii)  iRemnant == 1  then Mclu1 is a "soft" mass, whereas Mclu2 is a "normal" one;
  // iii) iRemnant == 2  then Mclu1 is a "normal" mass, whereas Mclu2 is a "soft" one;
  // iv)  otherwise           both Mclu1 and Mclu2 are "soft" masses.
  // Case i) is the standard one, always used for non-beam clusters;
  // case ii) and iii) occur for beam clusters with only one beam remnant and
  // when the flag _IOpRem is 1 (which is also its default value) which means 
  // that the cluster child containing the beam remnant should have "soft" mass
  // whereas the other one should have "normal" mass;
  // case iv), finally, occur for beam clusters with two beam remnants 
  // (regardless of the value of _IOpRem) or when _IOpRem is 0 which means 
  // that both cluster children masses should be soft.
  // Notice that, if the soft underlying event is on, then the beam clusters
  // are tagged as not available, therefore they will not be treated by
  // Cluster Fissioner (and the other Cluster Hadronization classes as well).
  // Finally, sometimes, when the phase space available is tiny, many attempts 
  // fail to produce a pair of masses kinematically acceptable, and in these cases 
  // it gives up returning false; otherwise it returns true when the splitting succeeds.

  void calculateKinematics(const Lorentz5Momentum & pClu, const Lorentz5Momentum & p0Q1, 
			   const bool decayOneHadronClu1, const bool decayOneHadronClu2,
			   Lorentz5Momentum & pClu1, Lorentz5Momentum & pClu2,      
			   Lorentz5Momentum & pQ1, Lorentz5Momentum & pQbar,      
			   Lorentz5Momentum & pQ, Lorentz5Momentum & pQ2bar ) const;
  // Determine the full kinematics of the fission of an heavy cluster C -> C1 + C2
  
  void calculatePositions( const Lorentz5Momentum & pClu, const LorentzPoint & positionClu,
			   const Lorentz5Momentum & pClu1, const Lorentz5Momentum & pClu2, 
			   LorentzPoint & positionClu1, LorentzPoint & positionClu2 ) const;
  // Determine the positions of the two children clusters.

  HadronsSelectorPtr _hadronsSelector;
  GlobParamPtr       _globalParameters;  

  Energy _ClMax;
  double _ClPow;
  double _PSplt1;
  double _PSplt2;

  Energy _BtClM; // At the moment it is not an interfaced parameter.
  int _IOpRem;   // At the moment it is not an interfaced parameter.

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of ClusterFissioner.
template <>
struct BaseClassTrait<Herwig::ClusterFissioner,1> {
  typedef Pythia7::HandlerBase NthBase;
};

// The following template specialization informs Pythia7 about the
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

#include "ClusterFissioner.icc"

#endif /* HERWIG_ClusterFissioner_H */
