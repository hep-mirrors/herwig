// -*- C++ -*-
#ifndef HERWIG_ClusterHadronizationHandler_H
#define HERWIG_ClusterHadronizationHandler_H
//
// This is the declaration of the <!id>ClusterHadronizationHandler<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is the main driver of the Cluster Hadronization: it is 
// responsible for the proper handling of all other specific collaborating
// classes (<!class>PartonSplitter<!!class>, <!class>ClusterFinder<!!class>, 
// <!class>ColourReconnector<!!class>, 
// <!class>ClusterFissioner<!!class>, <!class>LightClusterDecayer<!!class>, 
// <!class>ClusterDecayer<!!class>) 
// and for the storing of the produced particles in the Event record.
// 
// Important implementation detail: the private member <!id>_collecCluPtr<!!id>
// is the collection of all Cluster class objects that are created during the 
// whole cluster hadronization. This collection is initially cleaned, and 
// then it is passed to the collaborating classes which are responsible to 
// fill it properly. <BR>
// The elements of this collection, which are cluster objects, are properly
// interconnected in such a way to store the full, complete information about
// their origin, evolution, and end. This allows both to get at any time any
// information we could need (even for future, extensive changes) and to fully
// debug the cluster hadronization. <BR>
//
// Notice that the access to the <!class>GlobalParameters<!!class> class 
// instance is provided only to allow non-interfaced and non-persistent classes
// (<!class>Cluster<!!class>) to access the global parameters and/or to drawn 
// random numbers. This is done in the run initialization, <!id>doinitrun()
// <!!id>, by setting static pointers defined in those non-interfaced and 
// non-persistent classes.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:GlobalParameters.html">GlobalParameters.h</a>, <BR>
// <a href="http:PartonSplitter.html">PartonSplitter.h</a>, <BR>
// <a href="http:ClusterFinder.html">ClusterFinder.h</a>, <BR>
// <a href="http:ColourReconnector.html">ColourReconnector.h</a>, <BR>
// <a href="http:ClusterFissioner.html">ClusterFissioner.h</a>, <BR>
// <a href="http:LightClusterDecayer.html">LightClusterDecayer.h</a>, <BR>
// <a href="http:ClusterDecayer.html">ClusterDecayer.h</a>, <BR>
// <a href="http:Cluster.html">Cluster.h</a>.
// 

#include "Pythia7/Handlers/HadronizationHandler.h"
#include "Herwig++/Config/GlobalParameters.h"
#include "PartonSplitter.h"
#include "ClusterFinder.h"
#include "ColourReconnector.h"
#include "ClusterFissioner.h"
#include "LightClusterDecayer.h"
#include "ClusterDecayer.h"
#include "Cluster.h"


namespace Herwig {


using namespace Pythia7;


class ClusterHadronizationHandler: public HadronizationHandler {

public:

  inline ClusterHadronizationHandler();
  inline ClusterHadronizationHandler(const ClusterHadronizationHandler &);
  virtual ~ClusterHadronizationHandler();
  // Standard ctors and dtor.

public:

  virtual void handle(PartialCollisionHandler & ch, const tPVector & tagged,
		      const Hint & hint) throw(Veto, Stop, Exception);
  // The main method which manages the all cluster hadronization.

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
  virtual void doinitrun(); 
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<ClusterHadronizationHandler> initClusterHadronizationHandler;
  // Describe a concrete class with persistent data.

  ClusterHadronizationHandler & operator=(const ClusterHadronizationHandler &);
  // Private and non-existent assignment operator.

  void printStep(tStepPtr ptrStep, const string & title);
  // Print the step for debugging.

  void debuggingInfo(PartialCollisionHandler & ch, ClusterVector &);
  // Print information about the final, complete collections of clusters
  // for debugging.
 
  GlobParamPtr           _globalParameters;
  PartonSplitterPtr      _partonSplitter;
  ClusterFinderPtr       _clusterFinder;
  ColourReconnectorPtr   _colourReconnector;
  ClusterFissionerPtr    _clusterFissioner;
  LightClusterDecayerPtr _lightClusterDecayer;
  ClusterDecayerPtr      _clusterDecayer; 
};


}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of ClusterHadronizationHandler.
template <>
struct BaseClassTrait<Herwig::ClusterHadronizationHandler,1> {
  typedef HadronizationHandler NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ClusterHadronizationHandler>: public ClassTraitsBase<Herwig::ClusterHadronizationHandler> {
  static string className() { return "/Herwig++/ClusterHadronizationHandler"; }
  // Return the class name.
  static string library() { return "libHwHadronization.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "ClusterHadronizationHandler.icc"

#endif /* HERWIG_ClusterHadronizationHandler_H */
