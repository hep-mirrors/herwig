// -*- C++ -*-
#ifndef HERWIG_Cluster_H
#define HERWIG_Cluster_H
//
// This is the declaration of the <!id>Cluster<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class represents a cluster, which is a colour singlet made usually <BR>
// of two components (quark-antiquark, quark-diquark, antiquark-antidiquark) <BR>
// or rarely by three components (quark-quark-quark, antiquark-antiquark-antiquark). <BR>
// A reference to the container with the pointers to its Components is provided. <BR>
// The class provides access to the pointer to: <BR>
// <UL>
//  <LI> the cluster father, in the case that the cluster it is a fission product <BR> 
//       of an heavy cluster; or, occasionally, in the case that the cluster is  <BR>
//       the "redefinition" of another cluster (his father), for example in the case <BR>
//       that its three quarks or anti-quarks components have been redefined as two <BR>
//       components (quark+diquark, or antiquark+antidiquark), or when the cluster <BR>
//       has been the partner of the momentum reshuffling necessary to conserve <BR>
//       energy-momentum when a light cluster decays into a single hadron. 
//  <LI> its children clusters (usually two), in the case the cluster is an heavy cluster <BR>
//       that undergoes to fission; or occasionally, in the case the cluster has been <BR>
//       "redefined" (re-interpreted) for example in the case that its three quarks or <BR>
//       anti-quarks components have been redefined as two components (quark+diquark, or <BR>
//       antiquark+antidiquark), or when its momentum has been reshuffled.  
//  <LI> its (eventual) reshuffling partner, necessary for energy-momentum <BR> 
//       conservation when light clusters are decayed into single hadron.
//  <LI> its Hadrons (Particle class objects) decays. <BR>
//       Notice that, differently from what is in Fortran Herwig, in the (rare) <BR> 
//       cases in which an heavy cluster undergoes to fission and one of its <BR>
//       two fission products is light (that is below the threshold to produce <BR>
//       the lightest pair of hadrons with proper flavours), such light <BR>
//       product is identified with the lightest single hadron (with proper flavours) <BR>
//       without create an cluster. Therefore, in these cases, the father (heavy) cluster <BR>
//       has one single cluster child and one single hadron child (rather than to have <BR>
//       two cluster children as usual).
// </UL> 
// The rationale behind the design of this class has been to allow to save <BR>
// completely all of the information on any transformation process which <BR>
// affect the clusters, from their formation through their entire life. <BR>
// This, in its turn, is motivated to allow any kind of future extension <BR>
// in the Cluster Hadronization model, to avoid the typical, albeit familiar, <BR>
// bad situation in which one does not have (because not stored somewhere!) <BR>
// the information that he would need... Furthermore, this pletora of information <BR>
// is quite useful as debugging tool, providing many possible consistent checks. <BR>
//
// Notice that in order to determine the cluster position from the positions <BR>
// of the components, <!id>Cluster<!!id> class needs some global parameters. <BR>
// Because <!id>Cluster<!!id> class is neither interfaced nor persistent, <BR>
// a static pointer to the <!class>GlobalParameters<!!class> class instance, <BR>
// where the global parameters are, is used. Such static pointer is <BR>
// set via the method <!id>setPointerGlobalParameters<!!id>, during the <BR>
// run initialization, <!id>doinitrun()<!!id>, in the steering class
// <!class>ClusterHadronizationHandler<!!class>. <BR>
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
#include "Component.h"
#include "Herwig++/Config/GlobalParameters.h"


namespace Herwig {


  using namespace Pythia7;
 
  class Cluster : public ReferenceCounted {

  public:

    Cluster(tPPtr part1, tPPtr part2, tPPtr part3 = tPPtr() );    
    Cluster(tPPtr part1, const long id2);
    Cluster(const long id1, const long id2);

    static inline void setPointerGlobalParameters( Ptr<GlobalParameters>::pointer 
						   inputPointerGlobalParameters );
    // Set the static pointer to the GlobalParameters object.
    // It is set in <!id>ClusterHadronizationHandler::doinitrun()<!!id>. 

    inline Cluster(const Cluster & x);
    inline ~Cluster();
    // Constructors and destructor.

    inline int numComponents() const;
    // Number of quark (diquark) constituents (normally two).    

    Energy sumConstituentMasses() const;
    // Sum of the constituent masses of the components of the cluster.    

    inline Energy mass() const;
    inline void mass(const Energy inputMass);
    // Access/set cluster mass.

    inline const Lorentz5Momentum & momentum() const;
    inline void momentum(const Lorentz5Momentum & momentum);
    // Cluster 5-momentum.

    inline const LorentzPoint & position() const;
    inline void position(const LorentzPoint & position);
    // Cluster 4-position.
    
  public:

    inline const CollecCompPtr & components() const;
    // Return (a const reference to) the collection of components.

    inline tCluPtr parentCluster() const;
    inline void parentCluster(const tCluPtr inputCluPtr); 
    // Access/Set the pointer to the parent cluster.

    inline tCluPtr reshufflingPartnerCluster() const;
    inline void reshufflingPartnerCluster(const tCluPtr inputCluPtr);
    // Access/set the pointer to the reshuffling partner cluster.

    inline const CollecCluPtr & childrenClusters() const;
    // Return (a const reference to) the collection of the children clusters.

    void addChildrenClusters(const tCluPtr child1, const tCluPtr child2 = tCluPtr() ); 
    // Add children clusters.

    inline const PVector & childrenHadrons() const;
    // Return (a const reference to the collection of the children hadrons.

    void addChildrenHadrons(const tPPtr had1, 
			    const tPPtr had2 = tPPtr(), const tPPtr had3 = tPPtr() );
    // Add children hadrons.  

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

    inline bool isStatusFinal() const;
    // Return true if the cluster has hadron children.

  private:

    Cluster & operator=(const Cluster &);
    // Private and non-existent assignment operator.

    void calculateP();
    // Calculate the 5-momentum vector of the cluster
    void calculateX();
    // Calculate the 4-position vector of the cluster

    static Ptr<GlobalParameters>::pointer _pointerGlobalParameters;

    bool              _isAvailable;
    Lorentz5Momentum  _momentum;
    LorentzPoint      _position;
    CollecCompPtr     _collecCompPtr;

    tCluPtr      _parentCluPtr;             // it must be transient to avoid cycles
    tCluPtr      _reshufflingPartnerCluPtr; // it must be transient to avoid cycles
    CollecCluPtr _collecChildCluPtr; 
    PVector      _collecChildHadPtr;
  
  };


} // end namespace Herwig


#include "Cluster.icc"

#endif // HERWIG_Cluster_H 
