// -*- C++ -*-
#ifndef HERWIG_ClusteringGuide_H
#define HERWIG_ClusteringGuide_H
//
// This is the declaration of the ClusteringGuide class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ClusteringGuide.fh"

#include "Clusterer.h"
#include "CKKWHardProcess.h"
#include "Partitioner.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * ClusteringGuide is used to guide a cascade history
 * reconstruction such that clusterings which will yield
 * a non-appropriate hard process are not considered
 * during cascade history reconstruction.
 *
 *@author Simon Plaetzer
 *
 * @see \ref ClusteringGuideInterfaces "The interfaces"
 * defined for ClusteringGuide.
 */
class ClusteringGuide: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ClusteringGuide();

  /**
   * Construct giving clusterers, hard processes,
   * a partitioner and a set of particles and the
   * parent of the node.
   */
  inline ClusteringGuide (const vector<tClustererPtr>&,
			  const vector<tCKKWHardProcessPtr>&,
			  Partitioner *,
			  const vector<ClusteringParticleData>&,
			  tClusteringGuidePtr);

  /**
   * The destructor.
   */
  virtual ~ClusteringGuide();
  //@}

  /**@name Access the members */
  //@{

  /**
   * Return the clusterers used.
   */
  inline vector<tClustererPtr> clusterers() const;

  /**
   * Return the hard processes used.
   */
  inline vector<tCKKWHardProcessPtr> processes() const;

  /**
   * Return the particles.
   */
  inline vector<ClusteringParticleData> particles() const;

  /**
   * Return the possible clusterings and subsequent guides
   * starting from this configuration.
   */
  inline map<ClusteringConfigurationPtr,ClusteringGuidePtr> clusterings () const;

  /**
   * Return the parent.
   */
  inline tClusteringGuidePtr parent () const;

  //@}

public:

  /**@name Methods used during clustering. */
  //@{

  /**
   * Return true, if no more clusterings are left.
   */
  inline bool reachedHard () const;

  /**
   * Return a pointer to the hard proccess, if reachedHard()
   */
  inline tCKKWHardProcessPtr hard () const;

  /**
   * Set a pointer to the hard proccess
   */
  inline void hard (tCKKWHardProcessPtr);

  //@}

public:

  /**
   * Generate the guide starting from the current set 
   * of particles. Return true, if a valid hard process
   * has been reached.
   */
  bool generateGuide ();

public:

  /**
   * Set the event generator
   */
  inline void eventGenerator(tEGPtr);

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

  void debugDump (ostream&, const string& feed="");

  void debugDumpNode (ostream&,const string& feed="",bool showc=false);

#endif

private:

  /**
   * The clusterers to be used
   */
  vector<tClustererPtr> _clusterers;

  /**
   * The hard processes to look for
   */
  vector<tCKKWHardProcessPtr> _hardProcesses;

  /**
   * A pointer to the Partitioner to be used.
   * Note: The Partitioner is owned by CascadeReconstructor.
   */
  Partitioner * _partitioner;

  /**
   * The current particle configuration.
   */
  vector<ClusteringParticleData> _particles;

  /**
   * The possible clusterings on the current
   * particle configuration mapped to ClusteringGuides
   * starting from the emerging configuration.
   */
  map<ClusteringConfigurationPtr,ClusteringGuidePtr> _clusterings;

  /**
   * In case no more clusterings are left, this is
   * the hard subprocess.
   */
  vector<ClusteringParticleData> _hardProcess;

  /**
   * The hard process object which terminated this guide.
   */
  tCKKWHardProcessPtr _hard;

  /**
   * The parent of this node
   */
  tClusteringGuidePtr _parent;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<ClusteringGuide> initClusteringGuide;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ClusteringGuide & operator=(const ClusteringGuide &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ClusteringGuide. */
template <>
struct BaseClassTrait<Herwig::ClusteringGuide,1> {
  /** Typedef of the first base class of ClusteringGuide. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ClusteringGuide class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ClusteringGuide>
  : public ClassTraitsBase<Herwig::ClusteringGuide> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ClusteringGuide"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ClusteringGuide is implemented. It may also include several, space-separated,
   * libraries if the class ClusteringGuide depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "ClusteringGuide.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ClusteringGuide.tcc"
#endif

#endif /* HERWIG_ClusteringGuide_H */
