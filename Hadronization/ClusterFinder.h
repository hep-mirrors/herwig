// -*- C++ -*-
#ifndef HERWIG_ClusterFinder_H
#define HERWIG_ClusterFinder_H
/*! \class Herwig::ClusterFinder ClusterFinder.h "Herwig++\Hadronization\ClusterFinder.h
 * \brief This class forms clusters from the partons produced in the Shower.
 * \author Philip Stephens
 * \author Alberto Ribon
 * \ingroup Hadronization
 *
 * This class scans through the particles in the event and produces a 
 * collection of clusters, defined as a colour-singlet combinations of 
 * colour-connected particles. There are no assumptions about the type 
 * (i.e. quark or diquark) or number of the component particles of the 
 * cluster; however, most of the time clusters are formed by quark-antiquark 
 * pairs. In special situations, such as baryon-violating processes in
 * R-nonconserved Susy, three quarks (or three antiquarks) could form a 
 * cluster. Because at the moment we don't know how to handle 3-component 
 * clusters (i.e. how to fission heavy ones, or how to decay clusters), we 
 * provide also a separate method, reduceToTwoComponents, which 
 * does the job of redefining these 3-component clusters as "normal" 
 * 2-component ones, simply by randomly considering two (anti-) quarks as a 
 * (anti-) diquark. Notice that if in the future the method 
 * reduceToTwoComponents is modified or even eliminated, the 
 * main method for finding clusters, formClusters, will not need 
 * any change.
 */

#include <ThePEG/Handlers/HandlerBase.h>
#include "CluHadConfig.h"


namespace Herwig {


using namespace ThePEG;

class ThePEG::PartialCollisionHandler;  // forward declaration

class ClusterFinder: public ThePEG::HandlerBase {

public:

  inline ClusterFinder();
  inline ClusterFinder(const ClusterFinder &);
  virtual ~ClusterFinder();

  void formClusters(tCollPtr collisionPtr, const StepPtr & pstep, 
		    ClusterVector & clusters) throw(Veto, Stop, Exception);
  /*!< This routine forms the clusters of the event.
   *
   * Form clusters starting from the list of particles in the event.
   * It also checks if the cluster is a beam cluster, that is if
   * at least one of its components is a beam remnant.
   */

  void reduceToTwoComponents(const StepPtr &, ClusterVector&) 
    throw(Veto, Stop, Exception);
  /*!< Reduces three component clusters into two components.
   *
   * For the eventual clusters that have three components 
   * (quark, quark, quark) or (antiquark, antiquark, antiquark),
   * it redefines them as "normal" clusters with two components:
   * (quark,diquark) or (antiquark,antidiquark), by a random drawing.
   * This could be eliminated or changed in the future.
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
  //!< Called during an update phase.
  inline virtual void doinit() throw(InitException);
  //!< Called during the initialization phase.
  inline virtual void dofinish();
  //!< Called at finish phase.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  //!< Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  //!< Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<ClusterFinder> initClusterFinder;
  //!< Describe a concrete class with persistent data.

  ClusterFinder & operator=(const ClusterFinder &);
  //!< Private and non-existent assignment operator.
};

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of ClusterFinder.
template <>
struct BaseClassTrait<Herwig::ClusterFinder,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ClusterFinder>: public ClassTraitsBase<Herwig::ClusterFinder> {
  static string className() { return "/Herwig++/ClusterFinder"; }
  // Return the class name.
  static string library() { return "libHwHadronization.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#endif // DOXYGEN

#include "ClusterFinder.icc"

#endif /* HERWIG_ClusterFinder_H */
