// -*- C++ -*-
#ifndef HERWIG_ClusterFinder_H
#define HERWIG_ClusterFinder_H
//
// This is the declaration of the <!id>ClusterFinder<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class scans through the particles in the event and produces a 
// collection of clusters, defined as a colour-singlet combinations of 
// colour-connected particles. There are no assumptions about the type 
// (i.e. quark or diquark) or number of the component particles of the 
// cluster (however most of the time clusters are formed by quark-antiquark 
// pairs; but in special situations, as baryon-violating processes in
// R-nonconserved Susy, three quarks (or three antiquarks) could form a 
// cluster). Because at the moment we don't know how to handle 3-component 
// clusters (i.e. how to fission heavy ones, or how to decay clusters), we 
// provide also a separate method, <!id>reduceToTwoComponents<!!id>, which 
// does the job of redefining these 3-component clusters as "normal" 
// 2-component ones, simply by randomly considering two (anti-) quarks as a 
// (anti-) diquark. Notice that if in the future the method 
// <!id>reduceToTwoComponents<!!id> is modified or even eliminated, the 
// main method for finding clusters, <!id>formClusters<!!id>, will not need 
// any change.
//

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
  // Standard ctors and dtor.

  void formClusters(tCollPtr collisionPtr, const StepPtr & pstep, 
		    ClusterVector & clusters) throw(Veto, Stop, Exception);
  // Form clusters starting from the list of particles in the event.
  // It also checks if the cluster is a beam cluster, that is if
  // at least one of its components is a beam remnant.

  void reduceToTwoComponents(const StepPtr &, ClusterVector&) 
    throw(Veto, Stop, Exception);
  // For the eventual clusters that have three components 
  // (quark, quark, quark) or (antiquark, antiquark, antiquark),
  // it redefines them as "normal" clusters with two components:
  // (quark,diquark) or (antiquark,antidiquark), by a random drawing.
  // This could be eliminated or changed in the future.

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

  static ClassDescription<ClusterFinder> initClusterFinder;
  // Describe a concrete class with persistent data.

  ClusterFinder & operator=(const ClusterFinder &);
  //  Private and non-existent assignment operator.
};

}

// CLASSDOC OFF

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

#include "ClusterFinder.icc"

#endif /* HERWIG_ClusterFinder_H */
