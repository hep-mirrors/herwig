// -*- C++ -*-
#ifndef HERWIG_Cluster_H
#define HERWIG_Cluster_H
//
// This is the declaration of the <!id>Cluster<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class represents a cluster, which is a colour singlet made usually
// of two components (quark-antiquark, quark-diquark, antiquark-antidiquark)
// or rarely by three components (quark-quark-quark, antiquark-antiquark-
// antiquark).
// A reference to the container with the pointers to its Components is 
// provided. <BR>
// The class provides access to the pointer to:
// <UL>
//  <LI> the cluster father, in the case that the cluster it is a fission 
//       product 
//       of an heavy cluster; or, occasionally, in the case that the cluster is
//       the "redefinition" of another cluster (his father), for example in 
//       the case
//       that its three quarks or anti-quarks components have been redefined 
//       as two
//       components (quark+diquark, or antiquark+antidiquark), or when the 
//       cluster
//       has been the partner of the momentum reshuffling necessary to conserve
//       energy-momentum when a light cluster decays into a single hadron. 
//  <LI> its children clusters (usually two), in the case the cluster is an 
//       heavy cluster 
//       that undergoes to fission; or occasionally, in the case the cluster 
//       has been "redefined" (re-interpreted) for example in the case that 
//       its three quarks or
//       anti-quarks components have been redefined as two components (quark+
//       diquark, or 
//       antiquark+antidiquark), or when its momentum has been reshuffled.  
//  <LI> its (eventual) reshuffling partner, necessary for energy-momentum 
//       conservation when light clusters are decayed into single hadron.
//  <LI> its Hadrons (Particle class objects) decays.
//       Notice that, differently from what is in Fortran Herwig, in the (rare)
//       cases in which an heavy cluster undergoes to fission and one of its
//       two fission products is light (that is below the threshold to produce
//       the lightest pair of hadrons with proper flavours), such light
//       product is identified with the lightest single hadron (with proper 
//       flavours) without create an cluster. Therefore, in these cases, the 
//       father (heavy) cluster has one single cluster child and one single 
//       hadron child (rather than to have two cluster children as usual).
// </UL> 
// The rationale behind the design of this class has been to allow to save 
// completely all of the information on any transformation process which 
// affect the clusters, from their formation through their entire life. 
// This, in its turn, is motivated to allow any kind of future extension 
// in the Cluster Hadronization model, to avoid the typical, albeit familiar,
// bad situation in which one does not have (because not stored somewhere!) 
// the information that he would need... Furthermore, this pletora of 
// information is quite useful as debugging tool, providing many possible
// consistent checks. <BR>
//
// Notice that in order to determine the cluster position from the positions
// of the components, <!id>Cluster<!!id> class needs some global parameters.
// Because <!id>Cluster<!!id> class is neither interfaced nor persistent, 
// a static pointer to the <!class>GlobalParameters<!!class> class instance, 
// where the global parameters are, is used. Such static pointer is 
// set via the method <!id>setPointerGlobalParameters<!!id>, during the
// run initialization, <!id>doinitrun()<!!id>, in the steering class
// <!class>ClusterHadronizationHandler<!!class>. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:Component.html">Component.h</a>, <BR>
// <a href="http:GlobalParameters.html">GlobalParameters.h</a>, <BR>
// <a hraf="http:ClusterHadronizationHandler.html">ClusterHadronizationHandler.h</a>.
// 

#include "CluHadConfig.h"
#include "Pythia7/Pointer/Ptr.h"
#include "Pythia7/Pointer/ReferenceCounted.h"
#include "Pythia7/Pointer/PtrTraits.h"
#include "Pythia7/Pointer/RCPtr.h"
#include "Herwig++/Config/GlobalParameters.h"
#include "Pythia7/EventRecord/Particle.h"
#include "Herwig++/Utilities/EnumParticles.h"
#include <iostream>

namespace Herwig {

  using namespace Pythia7;
 
  class Cluster : public Particle {

  public:

    Cluster();
    Cluster(tcEventPDPtr);
    Cluster(tPPtr part1, tPPtr part2, tPPtr part3 = tPPtr());    
    Cluster(const Cluster &);
    Cluster(const Particle &);
    virtual ~Cluster();

    inline void * operator new(size_t);
    inline void operator delete(void *, size_t);  

    static inline void setPointerGlobalParameters(GlobParamPtr gp);
    // Set the static pointer to the GlobalParameters object.
    // It is set in <!id>ClusterHadronizationHandler::doinitrun()<!!id>. 

    inline int numComponents() const;
    // Number of quark (diquark) constituents (normally two).    

    Energy sumConstituentMasses() const;
    // Sum of the constituent masses of the components of the cluster.    

    inline void setMass(const Energy inputMass);
    // Set cluster mass.

    tPPtr particle(int i) const;

    bool isPerturbative(int) const;
    void setPerturbative(int,bool);

    bool isBeamRemnant(int) const;
    void setBeamRemnant(int,bool);

    int clusterId() const { return _id; }
  public:

    inline tClusterPtr reshufflingPartnerCluster() const;
    inline void reshufflingPartnerCluster(const tClusterPtr inputCluPtr);
    // Access/set the pointer to the reshuffling partner cluster.

    bool isBeamCluster() const;
    void isBeamCluster(tPPtr part);
    // The first method return true if the cluster is a beam cluster, 
    // defined as a cluster with at least one beam remnant between its 
    // components; false otherwise. The second method set the component 
    // (if any) that points to "part" as a beam remnant.

    bool isAvailable() const;
    void isAvailable(bool inputAvailable);
    // Check/set whether the cluster can be processed or should be skipped,
    // the default being true. It is useful for ignoring the cluster beams 
    // when the cluster hadronization is not supposed to treat them, i.e. 
    // soft underlying event is on and it will take care of them. 

    inline bool isStatusInitial() const;
    // Return true if the cluster does not have cluster parent.

    inline bool isReadyToDecay() const;
    // Return true if the cluster does not have cluster children and 
    // it is not already decayed (i.e. it does not have hadron children)
    // (to be used only after the fission of heavy clusters).    

    inline bool isRedefined() const;
    // Return true if the cluster has one and only one cluster children
    // and no hadron children: that means either that its three quarks or 
    // anti-quarks components have been redefined as two components 
    // (quark+diquark, or antiquark+antidiquark), or that the cluster 
    // has been used as a partner for the momentum reshuffling necessary 
    // to conserve energy-momentum when a light cluster is decayed into
    // a single hadron (notice that this latter light cluster has 
    //  isRedefined()  false, because it has an hadron child).
    // In both cases, the unique cluster children is the new redefined 
    // cluster. The two cases can be distinguish by the next method.  

    inline bool hasBeenReshuffled() const;
    // Return true when it has a reshuffling partner.
    // Notice that a cluster can have  hasBeenreshuffled()  true but
    //  isRedefined()  false: this is the case of a light cluster
    // that decays into a single hadron.

    bool isStatusFinal() const;
    // Return true if the cluster has hadron children.

  private:
    virtual PPtr clone() const;
    virtual PPtr fullclone() const;
    //virtual void rebind(const EventTranslationMap &);

    Cluster & operator=(const Cluster &);
    // Private and non-existent assignment operator.

    void calculateP();
    // Calculate the 5-momentum vector of the cluster
    void calculateX();
    // Calculate the 4-position vector of the cluster

    bool initPerturbative(tPPtr p);
    // Determines whether constituent p is perturbative or not

    static GlobParamPtr _globalParameters;
    // This is needed to know whether a cluster is from a perturbative
    // quark or not

    static ClassDescription<Cluster> initCluster; 

    bool        _isAvailable;
    tClusterPtr _reshufflingPartner; // transient to avoid cycles
    ParticleVector _component;
    //   vector<Particle> _component;
    vector<bool> _isBeamRemnant;
    vector<bool> _isPerturbative;
    int _numComp;
    long _id;
  };

} // end namespace Herwig  

namespace Pythia7 {

template <>
struct BaseClassTrait<Herwig::Cluster,1> {
  typedef EventRecordBase NthBase;
};
 
template <>
struct ClassTraits<Herwig::Cluster>:
    public ClassTraitsBase<Herwig::Cluster> {
  static string className() { return "/Herwig++/Cluster"; }
  //  static TPtr create() {
  //  return TPtr::Create(Herwig::Cluster(Particle(tcEventPDPtr())));
  //}
};

}
#include "Cluster.icc"

#endif // HERWIG_Cluster_H 
