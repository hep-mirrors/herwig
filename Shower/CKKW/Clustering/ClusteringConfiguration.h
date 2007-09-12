// -*- C++ -*-
#ifndef HERWIG_ClusteringConfiguration_H
#define HERWIG_ClusteringConfiguration_H
//
// This is the declaration of the ClusteringConfiguration class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ClusteringConfiguration.fh"

#include "ClusteringParticle.h"
#include "Clusterer.fh"

#include "ClusteringGuide.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * ClusteringConfiguration stores data associated with
 * a possible clustering during cascade reconstruction.
 *
 *@author Simon Plaetzer
 *
 * @see \ref ClusteringConfigurationInterfaces "The interfaces"
 * defined for ClusteringConfiguration.
 */
class ClusteringConfiguration: public Interfaced {

  // The ClusteringGuide is a good friend
  friend class ClusteringGuide;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ClusteringConfiguration();

  /**
   * Construct giving particles to be clustered,
   * the emerging ones and the clusterer to perform
   * the clustering.
   */
  inline ClusteringConfiguration (const vector<ClusteringParticleData>&,
				  const vector<ClusteringParticleData>&,
				  ClusteringInteractionType::ClusteringInteractionType,
				  tClustererPtr);

  /**
   * The destructor.
   */
  virtual ~ClusteringConfiguration();
  //@}

public:

  /**@name Methods to access members */
  //@{

  /**
   * Return the index set for the particles to be
   * clustered.
   */
  inline vector<unsigned int> clusteringIndices () const;

  /**
   * Return the data of particles to be clustered
   */
  inline vector<ClusteringParticleData> toBeClustered () const;

  /**
   * Return the data of particles emerging
   * from the clustering.
   */
  inline vector<ClusteringParticleData> emergingFromClustering () const;

  /**
   * Return the interaction type of this clustering
   */
  inline int interaction () const;

  /**
   * Return the clusterer to perform this
   * clustering.
   */
  inline tClustererPtr clusterer () const;

  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

#ifdef HERWIG_DEBUG_CKKW

public:

  virtual void debugDump (ostream&);

#endif

private:

  /**
   * The indices of the particles to be clustered
   */
  vector<unsigned int> _clusteringIndices;

  /**
   * The precise configuration of particles to be clustered,
   * in the same ordering as resulting from _clusterIndices
   */
  vector<ClusteringParticleData> _toBeClustered;

  /**
   * The configuration of particles emerging from
   * the clustering
   */
  vector<ClusteringParticleData> _emergingFromClustering;

  /**
   * The interaction type this clustering is considered
   * to be done with.
   */
  int _interaction;

  /**
   * The clusterer to perform clusterings of this
   * type.
   */
  tClustererPtr _clusterer;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<ClusteringConfiguration> initClusteringConfiguration;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ClusteringConfiguration & operator=(const ClusteringConfiguration &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ClusteringConfiguration. */
template <>
struct BaseClassTrait<Herwig::ClusteringConfiguration,1> {
  /** Typedef of the first base class of ClusteringConfiguration. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ClusteringConfiguration class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ClusteringConfiguration>
  : public ClassTraitsBase<Herwig::ClusteringConfiguration> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ClusteringConfiguration"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ClusteringConfiguration is implemented. It may also include several, space-separated,
   * libraries if the class ClusteringConfiguration depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "ClusteringConfiguration.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ClusteringConfiguration.tcc"
#endif

#endif /* HERWIG_ClusteringConfiguration_H */
