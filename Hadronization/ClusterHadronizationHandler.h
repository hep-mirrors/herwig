// -*- C++ -*-
#ifndef HERWIG_ClusterHadronizationHandler_H
#define HERWIG_ClusterHadronizationHandler_H
/*! \class Herwig::ClusterHadronizationHandler ClusterHadronizationHandler.h "Herwig++\Hadronization\ClusterHadronizationHandler.h"
 * \brief Class that controls the cluster hadronization algorithm.
 * \author Philip Stephens
 * \author Alberto Ribon
 * \ingroup Hadronization
 *
 * This class is the main driver of the cluster hadronization: it is 
 * responsible for the proper handling of all other specific collaborating
 * classes PartonSplitter, ClusterFinder, ColourReconnector, ClusterFissioner, 
 * LightClusterDecayer, ClusterDecayer; 
 * and for the storing of the produced particles in the Event record.
 * 
 * Notice that the access to the GlobalParameters class 
 * instance is provided only to allow non-interfaced and non-persistent classes
 * (Cluster) to access the global parameters and/or to draw 
 * random numbers. This is done in the run initialization, doinitrun()
 * by setting static pointers defined in those non-interfaced and 
 * non-persistent classes.
 *
 * See also: GlobalParameters.h, PartonSplitter.h, ClusterFinder.h,
 * ColourReconnector.h, ClusterFissioner.h, LightClusterDecayer.h,
 * ClusterDecayer.h, Cluster.h.
 */ 

#include <ThePEG/Handlers/HadronizationHandler.h>
#include "Herwig++/Utilities/GlobalParameters.h"
#include "PartonSplitter.h"
#include "ClusterFinder.h"
#include "ColourReconnector.h"
#include "ClusterFissioner.h"
#include "LightClusterDecayer.h"
#include "ClusterDecayer.h"
#include "Cluster.h"


namespace Herwig {


using namespace ThePEG;


class ClusterHadronizationHandler: public HadronizationHandler {

public:

  inline ClusterHadronizationHandler();
  inline ClusterHadronizationHandler(const ClusterHadronizationHandler &);
  virtual ~ClusterHadronizationHandler();

public:

  virtual void handle(PartialCollisionHandler & ch, const tPVector & tagged,
		      const Hint & hint) throw(Veto, Stop, Exception);
  /*!< The main method which manages the all cluster hadronization.
   *
   * This routine directs "traffic". It determines which function is called
   * and on which particles/clusters. This function also handles the 
   * situation of vetos on the hadronization.
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
  inline virtual void doinit() throw(InitException);
  virtual void doinitrun(); 
  inline virtual void dofinish();

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  //!< Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  //!< Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<ClusterHadronizationHandler> initClusterHadronizationHandler;
  //!< Describe a concrete class with persistent data.

  ClusterHadronizationHandler & operator=(const ClusterHadronizationHandler &);
  //!< Private and non-existent assignment operator.

  void printStep(tStepPtr ptrStep, const string & title);
  //!< Print the step for debugging.

  void debuggingInfo(PartialCollisionHandler & ch, ClusterVector &);
  //!< Print information about the final, complete collections of clusters.
 
  GlobParamPtr           _globalParameters;
  //!< This is a pointer to a Herwig::GlobalParameters object.
  PartonSplitterPtr      _partonSplitter;
  //!< This is a pointer to a Herwig::PartonSplitter object.
  ClusterFinderPtr       _clusterFinder;
  //!< This is a pointer to a Herwig::ClusterFinder object.
  ColourReconnectorPtr   _colourReconnector;
  //!< This is a pointer to a Herwig::ColourReconnector object.
  ClusterFissionerPtr    _clusterFissioner;
  //!< This is a pointer to a Herwig::ClusterFissioner object.
  LightClusterDecayerPtr _lightClusterDecayer;
  //!< This is a pointer to a Herwig::LightClusterDecayer object.
  ClusterDecayerPtr      _clusterDecayer; 
  //!< This is a pointer to a Herwig::ClusterDecayer object.
};


}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of ClusterHadronizationHandler.
template <>
struct BaseClassTrait<Herwig::ClusterHadronizationHandler,1> {
  typedef HadronizationHandler NthBase;
};

// The following template specialization informs ThePEG about the
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

#endif // DOXYGEN

#include "ClusterHadronizationHandler.icc"

#endif /* HERWIG_ClusterHadronizationHandler_H */
