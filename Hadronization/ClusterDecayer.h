// -*- C++ -*-
#ifndef HERWIG_ClusterDecayer_H
#define HERWIG_ClusterDecayer_H
//
// This is the declaration of the <!id>ClusterDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class decays the "normal" clusters (that is the ones not too heavy <BR> 
// to be splitted, and not too light to decay into one hadron). <BR>
// These clusters decay into two hadrons, but it is straightforward <BR>
// to extend to the case of three body decays (in the case that one is <BR>
// insterested to explore this possibilities...). 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:HadronSelector.html">HadronSelector.h</a>.
// 

#include <ThePEG/Handlers/HandlerBase.h>
#include <ThePEG/EventRecord/Step.h>
#include "CluHadConfig.h"
#include "HadronSelector.h"
#include "Herwig++/Utilities/GlobalParameters.h"


namespace Herwig {


using namespace ThePEG;

  //class Cluster;             // forward declaration
class ThePEG::Particle;   // forward declaration


class ClusterDecayer: public ThePEG::HandlerBase {

public:

  inline ClusterDecayer();
  inline ClusterDecayer(const ClusterDecayer &);
  virtual ~ClusterDecayer();
  // Standard ctors and dtor.

  void decay(const StepPtr&) 
    throw(Veto, Stop, Exception);
  // Decays all clusters (not already decayed into a single hadron) into hadrons. 

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

  static ClassDescription<ClusterDecayer> initClusterDecayer;
  // Describe a concrete class with persistent data.

  ClusterDecayer & operator=(const ClusterDecayer &);
  //  Private and non-existent assignment operator.

public:
  pair<PPtr,PPtr> decayIntoTwoHadrons(tClusterPtr ptr) 
    throw(Veto, Stop, Exception);
  // It decays the cluster into two hadrons. 

private:
  void calculatePositions( const Lorentz5Momentum &, const LorentzPoint &, 
			   const Lorentz5Momentum &, const Lorentz5Momentum &,
			   LorentzPoint &, LorentzPoint &) const;
  // It calculates the positions of the children hadrons by 
  // gaussian smearing, with width inversely proportional to 
  // the cluster mass. around the parent cluster position.

  Ptr<HadronSelector>::pointer _hadronsSelector;
  Ptr<GlobalParameters>::pointer _globalParameters;
  
  int _ClDir1;
  int _ClDir2;
  double _ClSmr1;
  double _ClSmr2;

};


}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of ClusterDecayer.
template <>
struct BaseClassTrait<Herwig::ClusterDecayer,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ClusterDecayer>: public ClassTraitsBase<Herwig::ClusterDecayer> {
  static string className() { return "/Herwig++/ClusterDecayer"; }
  // Return the class name.
  static string library() { return "libHwHadronization.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "ClusterDecayer.icc"

#endif /* HERWIG_ClusterDecayer_H */
