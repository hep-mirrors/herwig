// -*- C++ -*-
#ifndef HERWIG_EmitterSpectatorConfiguration_H
#define HERWIG_EmitterSpectatorConfiguration_H
//
// This is the declaration of the EmitterSpectatorConfiguration class.
//

#include "Herwig++/Shower/CKKW/Clustering/ClusteringConfiguration.h"
#include "EmitterSpectatorConfiguration.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * EmitterSpectatorConfiguration is a specific implementation
 * of ClusteringConfiguration for algorithms working in an
 * emitter-spectator scenario.
 *
 *@author Simon Plaetzer
 *
 * @see \ref EmitterSpectatorConfigurationInterfaces "The interfaces"
 * defined for EmitterSpectatorConfiguration.
 */
class EmitterSpectatorConfiguration: public ClusteringConfiguration {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline EmitterSpectatorConfiguration();

  /**
   * Construct giving particles to be clustered,
   * the emerging ones and the clusterer to perform
   * the clustering. In addition, give the spectator's
   * indices (see EmitterSpectatorClustering)
   */
  inline EmitterSpectatorConfiguration (const vector<ClusteringParticleData>&,
					const unsigned int,
					const vector<ClusteringParticleData>&,
					const unsigned int,
					ClusteringInteractionType::ClusteringInteractionType,
					tClustererPtr);

  /**
   * The destructor.
   */
  virtual ~EmitterSpectatorConfiguration();
  //@}

public:

  /**
   * Return the emitter
   */
  inline ClusteringParticleData emitter() const;

  /**
   * Return the emission products
   */
  inline pair<ClusteringParticleData,ClusteringParticleData> emission () const;

  /**
   * Return the spectator before clustering
   */
  inline ClusteringParticleData spectatorBeforeClustering () const;

  /**
   * Return the spectator after clustering
   */
  inline ClusteringParticleData spectatorAfterClustering () const;

  /**
   * Return the spectator's index after emission
   */
  inline unsigned int spectatorBeforeClusteringIndex () const;

  /**
   * Return the emitter's index before emission
   */
  inline unsigned int emitterAfterClusteringIndex () const;

  /**
   * Return the spectator's index before emission
   */
  inline unsigned int spectatorAfterClusteringIndex () const;

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
   * The index of the spectator in the childrens vector
   */
  unsigned int _childSpectator;

  /**
   * The index of the spectator in the parents vector
   */
  unsigned int _parentSpectator;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<EmitterSpectatorConfiguration> initEmitterSpectatorConfiguration;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EmitterSpectatorConfiguration & operator=(const EmitterSpectatorConfiguration &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of EmitterSpectatorConfiguration. */
template <>
struct BaseClassTrait<Herwig::EmitterSpectatorConfiguration,1> {
  /** Typedef of the first base class of EmitterSpectatorConfiguration. */
  typedef Herwig::ClusteringConfiguration NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EmitterSpectatorConfiguration class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::EmitterSpectatorConfiguration>
  : public ClassTraitsBase<Herwig::EmitterSpectatorConfiguration> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::EmitterSpectatorConfiguration"; }
  /**
   * The name of a file containing the dynamic library where the class
   * EmitterSpectatorConfiguration is implemented. It may also include several, space-separated,
   * libraries if the class EmitterSpectatorConfiguration depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "EmitterSpectatorConfiguration.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EmitterSpectatorConfiguration.tcc"
#endif

#endif /* HERWIG_EmitterSpectatorConfiguration_H */
