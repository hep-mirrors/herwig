// -*- C++ -*-
#ifndef HERWIG_Cluster_H
#define HERWIG_Cluster_H

#include <ThePEG/Pointer/Ptr.h>
#include <ThePEG/Pointer/PtrTraits.h>
#include <ThePEG/Pointer/RCPtr.h>
#include <ThePEG/EventRecord/Particle.h>
#include "Herwig++/Utilities/EnumParticles.h"
#include <iostream>
#include "CluHadConfig.h"
#include "ClusterHadronizationHandler.fh"
#include "Cluster.fh"

namespace Herwig {
using namespace ThePEG;
 
/** \ingroup Hadronization
 *  \class Cluster
 *  \brief This class describes a cluster object.
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *
 *  This class represents a cluster, which is a colour singlet made usually
 *  of two components (quark-antiquark, quark-diquark, antiquark-antidiquark)
 *  or rarely by three components (quark-quark-quark, antiquark-antiquark-
 *  antiquark). A reference to the container with the pointers to its 
 *  Components is provided.
 *
 *  The class provides access to the pointers which point to:
 *
 *     - The cluster parent. In the case that the cluster it is a fission 
 *       product of a heavy cluster the parent is a cluster. If the cluster
 *       is formed from the perturbative partons then the parents will be
 *       the colour connected partons that formed the cluster.
 *     - The children (usually two). In the case the cluster is a
 *       heavy cluster that undergoes fission the children are clusters.
 *       Occasionally the cluster has been "redefined" (re-interpreted). For 
 *       example in the case that three quark or anti-quark components 
 *       have been redefined as two components (quark+diquark, or antiquark+
 *       antidiquark).
 *     - The (eventual) reshuffling partner, necessary for energy-momentum 
 *       conservation when light clusters are decayed into single hadron. Not
 *       all clusters will have a reshuffling partner.
 *  
 *  Notice that in order to determine the cluster position from the positions
 *  of the components, the Cluster class needs some parameters.
 *  Because the Cluster class is neither interfaced nor persistent, 
 *  a static pointer to the ClusterHadronizationHandler class instance, 
 *  where the parameters are, is used. This static pointer is 
 *  set via the method setPointerClusterHadHandler(), during the
 *  run initialization, doinitrun() of ClusterHadronizationHandler.
 *
 *  @see ClusterHadronizationHandler
 */ 
class Cluster : public Particle {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  Cluster();
  
  /**
   * Constructor with a particleData pointer
   */
  Cluster(tcEventPDPtr);
  
  /**
   * This creates a cluster from 2 (or 3) partons.
   */
  Cluster(tPPtr part1, tPPtr part2, tPPtr part3 = tPPtr());
  
  /**
   * Also a constructor where a particle is given not a cluster.
   */
  Cluster(const Particle &);

  /**
   * Particle uses the FixedSizeAllocator for (de)allocation.
   */
  inline void * operator new(size_t);
  
  /**
   * Particle uses the FixedSizeAllocator for (de)allocation.
   */
  inline void operator delete(void *, size_t);
  //@}

  /**
   * Set the static pointer to the ClusterHadronizationHandler object.
   * The pointer is set in ClusterHadronizationHandler::doinitrun().
   */ 
  static inline void setPointerClusterHadHandler(tcCluHadHdlPtr gp);
  
  /**
   * Number of quark (diquark) constituents (normally two).    
   */
  inline int numComponents() const;
  
  /**
   * Sum of the constituent masses of the components of the cluster.    
   */
  Energy sumConstituentMasses() const;
  
  /**
   * Set cluster mass.
   */
  inline void setMass(const Energy inputMass);
  
  /**
   * Returns the ith constituent.
   */
  tPPtr particle(int i) const;
  
  /**
   * Returns whether the ith constituent is from a perturbative process.
   */
  bool isPerturbative(int) const;
  
  /**
   * Sets whether the ith constituent is from a perturbative process.
   */
  void setPerturbative(int,bool);
  
  /**
   * Indicates whether the ith constituent is a beam remnant.
   */
  bool isBeamRemnant(int) const;
  
  /**
   * Sets whether the ith constituent is a beam remnant.
   */
  void setBeamRemnant(int,bool);
  
  /**
   * Returns the clusters id, not the same as the PDG id.
   */
  int clusterId() const { return _id; }

public:

  /**
   * Returns the cluster that was used for reshuffling.
   */
  inline tClusterPtr reshufflingPartnerCluster() const;
  
  /**
   * Set the pointer to the reshuffling partner cluster.
   */
  inline void reshufflingPartnerCluster(const tClusterPtr inputCluPtr);
  
  /**
   * Returns true when a constituent is a beam remnant.
   */
  bool isBeamCluster() const;
  
  /**
   * Sets the component (if any) that points to "part" as a beam remnant.
   */
  void isBeamCluster(tPPtr part);
  
  /**
   * Returns true if this cluster is to be handled by the hadronization.
   */
  bool isAvailable() const;
  
  /**
   * Sets the value of availability. 
   */
  void isAvailable(bool inputAvailable);
  
  /**
   * Return true if the cluster does not have cluster parent.
   */
  inline bool isStatusInitial() const;
  
  /** 
   * Return true if the cluster does not have cluster children and 
   * it is not already decayed (i.e. it does not have hadron children)
   * (to be used only after the fission of heavy clusters).    
   */
  inline bool isReadyToDecay() const;
  
  /**
   * Return true if the cluster has one and only one cluster children
   * and no hadron children: that means either that its three quarks or 
   * anti-quarks components have been redefined as two components 
   * (quark+diquark, or antiquark+antidiquark), or that the cluster 
   * has been used as a partner for the momentum reshuffling necessary 
   * to conserve energy-momentum when a light cluster is decayed into
   * a single hadron (notice that this latter light cluster has 
   *  isRedefined()  false, because it has an hadron child).
   * In both cases, the unique cluster children is the new redefined 
   * cluster. The two cases can be distinguish by the next method.  
   */
  inline bool isRedefined() const;
  
  /**
   * Return true when it has a reshuffling partner.
   * Notice that a cluster can have  hasBeenReshuffled()  true but
   *  isRedefined()  false: this is the case of a light cluster
   * that decays into a single hadron.
   */
  inline bool hasBeenReshuffled() const;
  
  /**
   * Return true if the cluster has hadron children.
   */
  bool isStatusFinal() const;
  
protected:
  
  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual PPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual PPtr fullclone() const;
  //@}

private:   
  /**
   * Private and non-existent assignment operator.
   */
  Cluster & operator=(const Cluster &);
  
  /**
   * Calculate the 5-momentum vector of the cluster
   * The 5-momentum of the cluster is given by
   * \f[ P = \sum_i p_i \f]
   * and the mass of the cluster is \f$m^2 = P^2\f$
   */
  void calculateP();
  
  /**
   * Calculate the 4-position vector of the cluster
   * Displacement of the ith constituent given by momentum \f$p_i\f$
   * vertex \f$x_i\f$ and mass \f$m_i\f$ is 
   * \f[ D_i = -C \log(r) \frac{p_i}{\sqrt{(p_i^2 - m_i^2)^2 + v^4}} \f]
   * where \f$r\f$ is a random number [0,1], 
   * \f$v\f$ is the minimum virtuality and \f$C\f$ is 
   * a conversion factor from GeV to millimeters. We can then find the 
   * difference in \f$s\f$ factors as
   * \f[ (s_1-s_2) = \frac{(\vec{p}_1 + \vec{p}_2) \cdot (\vec{x}_2 - 
   *                 \vec{x}_1)}{(\vec{p}_1 + \vec{p}_2) \cdot \vec{D}_1}. 
   * \f]
   * if \f$s_2>s_1\f$ then \f$s_1 = 1.0\f$ otherwise \f$s_2 = 1.0\f$. 
   * These are then used to determine the value of the clusters vertex as
   * \f[ X = \frac{1}{2} ( x_1 +x_2 + s_1 D_1 + s_2 D_2). \f]
   */
  void calculateX();
  
  /**
   * Determines whether constituent p is perturbative or not.
   */
  bool initPerturbative(tPPtr p);
  
  /**
   * This is needed to determine if a cluster is from a perturbative quark.
   */
  static tcCluHadHdlPtr _clusterHadHandler;
  
  /**
   * Describe an abstract base class with persistent data.
   */
  static ClassDescription<Cluster> initCluster; 
  
  bool        _isAvailable;        //!< Whether the cluster is hadronizing
  tClusterPtr _reshufflingPartner; //!< The partner used in reshuffling
  ParticleVector _component;       //!< The constituent partons
  tParticleVector _original;       //!< The original components
  vector<bool> _isBeamRemnant;     //!< Whether a parton is a beam remnant
  vector<bool> _isPerturbative;    //!< Whether a parton is perturbative
  int _numComp;                    //!< The number of constituents
  long _id;                        //!< The id of this cluster
};
  
} // end namespace Herwig  

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of Cluster.
 */
template <>
struct BaseClassTrait<Herwig::Cluster,1> {
  /** Typedef of the base class of Cluster. */
  typedef Particle NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */ 
template <>
struct ClassTraits<Herwig::Cluster>:
  public ClassTraitsBase<Herwig::Cluster> {
  /** Return the class name. */
  static string className() { return "Herwig::Cluster"; }
};

/** @endcond */

}

#include "Cluster.icc"

#endif // HERWIG_Cluster_H 
