// -*- C++ -*-
#ifndef HERWIG_Component_H
#define HERWIG_Component_H
//
// This is the declaration of the <!id>Component<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// Objects of this class represent one component (constituent) <BR>
// of a cluster, that is a (anti-)quark or a (anti-)diquark.   <BR>
// In the case of clusters produced by heavy cluster fission, at least <BR>
// one of the cluster component does not correspond to a "Particle", <BR>
// but it is simply defined by its flavour (that can be inferred <BR>
// from the PDG id number corresponding to the quark or diquark <BR>
// which carries the same flavour). For the other "normal" components, <BR>
// the class <!id>Component<!!id> is a simple wrap around Particle, and it offers <BR>
// access to the Particle object via the method <I> pointerParticle()</I>.<BR>
// Notice that the momentum of the Component is kept as a (private) <BR>
// member of the class, and is set equal to the momentum of the <BR>
// pointed particle only in the constructor, but later, the two <BR>
// momenta need not to be the same. For example, in the class
// <!class>LightClusterDecayer<!!class> <BR>, 
// to preserve the full information before and after the reshuffling of momenta <BR>
// (necessary to preserve energy-momentum in the case of decay of light <BR>
// clusters into single hadrons) the momentum of the component is different <BR>
// with respect to the one of the pointed particle: the latter corresponds to the <BR>
// value before the reshuffling, whereas the former corresponds to the value after it. <BR>
// Another similar situation is in the class <!class>ClusterFissioner<!!class>. <BR>
// The mass of the component is defined as the 5-th component of the 5-momentum <BR>
// associated to it. Such mass can be incosistent (different) with respect to the <BR>
// invariant mass obtained from the 4 "normal" components of the momentum. <BR>
// This is not harmful, because it is proper handled in all the classes that <BR>
// deal with components. <BR>
//
// Notice: <BR> 
// <UL>
//  <LI> In the constructor <!id>Component(tPPtr pptr)<!!id> the <BR>
//       effective constituent gluon mass, which a global parameter defined in <BR>
//       <!class>GlobalParameters<!!class>, is needed. <BR>
//       Because <!id>Component<!!id> class is neither interfaced nor persistent, <BR>
//       a static pointer to the <!id>GlobalParameters<!!id> class instance is used. <BR>
//       Such static pointer is set via the method <!id>setPointerGlobalParameters<!!id>, <BR>
//       during the run initialization, <!id>doinitrun()<!!id>, in the steering class <BR>
//       <!class>ClusterHadronizationHandler<!!class>. 
// </UL>
// 
// CLASSDOC SUBSECTION See also:
//
// <a href="http:CluHadConfig.html">CluHadConfig.h</a>, <BR>
// <a href="http:GlobalParameters.html">GlobalParameters.h</a>, <BR>
// <a hraf="http:ClusterHadronizationHandler.html">ClusterHadronizationHandler.h</a>.
// 

#include "CluHadConfig.h"
#include "Pythia7/Pointer/Ptr.h"
#include "Pythia7/Pointer/ReferenceCounted.h"
#include "Pythia7/Pointer/PtrTraits.h"
#include "Pythia7/Pointer/RCPtr.h"
#include "Pythia7/EventRecord/Particle.h"
#include "Herwig++/Config/GlobalParameters.h"


namespace Herwig {
 
  
  using namespace Pythia7;

  class Component : public ReferenceCounted {

  public:

    explicit inline Component(const long id);
    explicit Component(tPPtr pptr);
    inline Component(const Component & x);
    inline ~Component();
    // Constructors and Destructor.

    static inline void setPointerGlobalParameters( Ptr<GlobalParameters>::pointer 
						   inputPointerGlobalParameters );
    // Set the static pointer to the GlobalParameters object.
    // It is set in <!id>ClusterHadronizationHandler::doinitrun()<!!id>.
						    
    inline long id() const;
    inline void id(const long inputIn);
    // Access/Set the id of the component (parton or diquark).

    inline const Energy mass() const;
    inline void mass(const Energy & inputMass);
    // Access/Set the mass of the component (parton or diquark).

    inline tPPtr pointerParticle() const;
    inline void pointerParticle(tPPtr inputPptr);
    // Access/Set the pointer to the Particle object associated to the 
    // component (parton or diquark).

    inline bool isPerturbative() const;
    inline void isPerturbative(const bool inputPerturbative);
    // Access/Set the logical that is true is the parton has perturbative 
    // origin (either from the hard subprocess or from showering).

    inline bool isBeamRemnant() const;
    inline void isBeamRemnant(const bool inputBeamRemnant);
    // Access/Set the logical that is true is the parton is a beam remnant.
    // The default is false. (This flag could be avoided by using 
    // pointerParticle(), but it is more efficient to use it).

    inline const Lorentz5Momentum & momentum() const;
    inline void momentum(const Lorentz5Momentum & inputMomentum);
    // Access/Set the 5-momentum vector of the component (parton or diquark).

    inline const LorentzPoint & position() const;
    inline void position(const LorentzPoint & inputPosition);
    // Return the 4-position vector of the component (parton or diquark).
    
  private:
    
    Component & operator=(const Component &);
    // Private and non-existent assignment operator.
    
    static Ptr<GlobalParameters>::pointer _pointerGlobalParameters;

    PPtr _pptr;
    long _id;
    bool _perturbative;
    bool _beamRemnant;
    Lorentz5Momentum _momentum;
    LorentzPoint _position;

  };
  
  
} // end namespace Herwig


#include "Component.icc"

#endif // HERWIG_Component_H 
