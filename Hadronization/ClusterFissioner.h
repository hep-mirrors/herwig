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
// children are draw from a power-like mass distribution; for beam clusters,
// according to the value of the flag <!id>_IOpRem<!!id>, either both 
// children masses are draw from a fast-decreasing exponential mass 
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
// <a href="http:HadronSelector.html">HadronSelector.h</a>.
// 

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
  // Split either heavy clusters or beam clusters recursively until all 
  // children have mass below some threshold. 
  // For beam clusters, they are split only if the soft underlying event
  // is switched off, otherwise these clusters will be tagged as unavailable
  // and they will not be treated by the hadronization altogether. 
  // In the case beam clusters will be split, the procedure is exactly
  // the same as for normal non-beam clusters, with the only exception
  // of the mass spectrum from which to draw the masses of the two 
  // cluster children (see method <!id>drawChildrenMasses<!!id> for details).
    
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

  void cut(tClusterPtr, const StepPtr&, ClusterVector&);
public:
  typedef pair<PPtr,PPtr> PPair;
  typedef pair<PPair,PPair> cutType;
  cutType cut(tClusterPtr &);
  // Split the input cluster (which can be either an heavy non-beam
  // cluster or a beam cluster). The result is two pairs of particles. The
  // first element of each pair is new cluster/hadron, while the second
  // element of each pair is the particle drawn from the vacuum to create
  // the new cluster/hadron.
  // Notice that this method treats also beam clusters by using a different
  // mass spectrum used to generate the cluster child masses (see method
  // drawChildMass).

private:
  PPair produceHadron(const long id1, const long id2, Lorentz5Momentum &a,
		      LorentzPoint &b) const;
  PPair produceCluster(tPPtr &p1, const long id, Lorentz5Momentum &a, 
		       LorentzPoint &b, Lorentz5Momentum &c, 
		       Lorentz5Momentum &d, const bool rem) const;
  // This routine produces a new cluster with the flavours given by p1 and id.
  // The new 5 momentum is a and the parent momentum are c and d. C is for the
  // p1 and d is for the new particle id. rem specifies whether the existing
  // particle is a beam remnant or not.

  long drawNewFlavour() const;
  // Return the id ( >0 ) of the quark-antiquark pair from the vacuum
  // needed for fission of an heavy cluster. Equal probabilities
  // are assumed for quarks  u , d , s . 

  void drawChildMass(const Energy M, const Energy m1, const Energy m2, 
		     const Energy m3, Energy & Mc, const double exp,
                           const double b, const bool rem) const; 
  // Draw the masses Mc of the the cluster child produced 
  // by the fission of an heavy cluster (of mass M). m1, m2 are the masses
  // of the constituents of the cluster; m3 is the mass of the quark extract 
  // from the vacuum (together with its antiparticle). The algorithm produces
  // the mass of the cluster formed with consituent m1.
  // Two mass distributions can be used for the child cluster mass:
  // 1) power-like mass distribution ("normal" mass) with power exp 
  // 2) fast-decreasing exponential mass distribution ("soft" mass) with 
  //    rmin. rmin is given by exp(-2*maxPS/average) where average is the
  //    given by the average method parameter.
  // The choice of which mass distribution should be used for each of the two
  // cluster children is dictated by the bool rem. If _IOpRem is 0, the
  // soft distribution is always used.
  // Finally, sometimes, when the phase space available is tiny, many attempts 
  // fail to produce a pair of masses kinematically acceptable; in these cases 
  // it gives up returning false, otherwise it returns true when the splitting 
  // succeeds.

  void calculateKinematics(const Lorentz5Momentum &pClu, 
		           const Lorentz5Momentum &p0Q1, 
			   const bool toHadron1, const bool toHadron2,
			   Lorentz5Momentum &pClu1, Lorentz5Momentum &pClu2, 
			   Lorentz5Momentum &pQ1, Lorentz5Momentum &pQb, 
			   Lorentz5Momentum &pQ2, Lorentz5Momentum &pQ2b) const;
  // Determine the full kinematics of the fission of an heavy cluster 
  // C -> C1 + C2
  
  void calculatePositions(const Lorentz5Momentum &pClu, 
		          const LorentzPoint & positionClu,
			  const Lorentz5Momentum & pClu1, 
			  const Lorentz5Momentum & pClu2, 
			  LorentzPoint & positionClu1, 
			  LorentzPoint & positionClu2 ) const;
  // Determine the positions of the two children clusters.

  HadronSelectorPtr _hadronsSelector;
  GlobParamPtr      _globalParameters;  
  Energy _ClMax;
  double _ClPow;
  double _PSplt1;
  double _PSplt2;

  Energy _BtClM; // At the moment it is not an interfaced parameter.
  int _IOpRem;   // At the moment it is not an interfaced parameter.

};

}

// CLASSDOC OFF

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

#include "ClusterFissioner.icc"

#endif /* HERWIG_ClusterFissioner_H */
