// -*- C++ -*-
#ifndef HERWIG_LightClusterDecayer_H
#define HERWIG_LightClusterDecayer_H
//
// This is the declaration of the <!id>LightClusterDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This is the class that performs the decay of light clusters into <BR>
// only one hadron. The major difficulty is that a kinematical reshuffling <BR>
// is necessary, between the cluster under consideration and its <BR>
// "neighbouring" clusters, to conserve energy-momentum in one-body decay. <BR>
// Notice that, differently from what happens in Fortran Herwig, <BR>
// light (that is below the threshold for the production of the lightest <BR>
// pair of hadrons with the proper flavours) fission products, produced <BR>
// by the fission of heavy clusters in <!class>ClusterFissioner<!!class> <BR>
// have been already "decayed" into single hadron (the lightest one <BR>
// with proper flavour) by the same latter class, without require <BR>
// any reshuffling. Therefore the light clusters that are treated in <BR>
// this <!id>LightClusterDecayer<!!id> class are produced directly <BR>
// (originally) by the <!class>ClusterFinder<!!class>. <BR>
//	
// Notice:
// <UL>
//  <LI> The choice of the candidate cluster with whom to reshuffle momentum <BR>
//       is based on the minimal space-time distance from the light cluster.
// </UL>
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:HadronsSelector.html">HadronsSelector.h</a>.
// 

#include "Pythia7/Handlers/HandlerBase.h"
#include "CluHadConfig.h"
#include "HadronsSelector.h"


namespace Herwig {


using namespace Pythia7;

class Cluster; // forward declaration


class LightClusterDecayer: public Pythia7::HandlerBase {

public:

  inline LightClusterDecayer();
  inline LightClusterDecayer(const LightClusterDecayer &);
  virtual ~LightClusterDecayer();
  // Standard ctors and dtor.

  void decay(CollecCluPtr & collecCluPtr) throw (Veto, Stop, Exception);
  // It does the decay of light hadron in one hadron: this requires
  // also a kinematical reshuffling for energy-momentum conservation
  // (this is done explicitly by the (private) method reshuffling() ).

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

  static ClassDescription<LightClusterDecayer> initLightClusterDecayer;
  // Describe a concrete class with persistent data.

  LightClusterDecayer & operator=(const LightClusterDecayer &);
  //  Private and non-existent assignment operator.

  bool reshuffling( const long idhad1, tCluPtr cluPtr1, tCluPtr cluPtr2 ,
		    vector<CluPtr> & vecNewRedefinedCluPtr ) 
    throw (Veto, Stop, Exception); 
  // This (private) method, called by decay(), takes care of the kinematical
  // reshuffling necessary for energy-momentum conservation.
  
  Ptr<HadronsSelector>::pointer _pointerHadronsSelector;

  double _B1Lim;

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of LightClusterDecayer.
template <>
struct BaseClassTrait<Herwig::LightClusterDecayer,1> {
  typedef Pythia7::HandlerBase NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::LightClusterDecayer>: public ClassTraitsBase<Herwig::LightClusterDecayer> {
  static string className() { return "/Herwig++/LightClusterDecayer"; }
  // Return the class name.
  static string library() { return "libHwHadronization.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "LightClusterDecayer.icc"

#endif /* HERWIG_LightClusterDecayer_H */
