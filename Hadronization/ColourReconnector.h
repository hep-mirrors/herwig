// -*- C++ -*-
#ifndef HERWIG_ColourReconnector_H
#define HERWIG_ColourReconnector_H
//
// This is the declaration of the <!id>ColourReconnector<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class does the nonperturbative colour rearrangement, after the <BR>
// nonperturbative gluon splitting and the "normal" cluster formation. <BR>
// It uses the list of particles in the event record, and the collections of <BR>
// "usual" clusters which is passed to the main method. If the colour reconnection <BR>
// is actually accepted, then the previous collections of "usual" <BR>
// clusters is first deleted and then the new one is created.
//
// NB) This class must be implemented!
//

#include "Pythia7/Handlers/HandlerBase.h"
#include "CluHadConfig.h"


namespace Herwig {


using namespace Pythia7;

class Pythia7::PartialCollisionHandler; // forward declaration


class ColourReconnector: public Pythia7::HandlerBase {

public:

  inline ColourReconnector();
  inline ColourReconnector(const ColourReconnector &);
  virtual ~ColourReconnector();
  // Standard ctors and dtor.

  void rearrange(PartialCollisionHandler & ch, const StepPtr & pstep,
                 ClusterVector & clusters) throw(Veto, Stop, Exception);
  // Does the colour rearrangement, starting from the list of particles
  // in the event record, and the collection of "usual" clusters passed
  // in input. If the actual rearrangement is accepted, the new collection 
  // of clusters is overriden to the intial one.
    
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

  static ClassDescription<ColourReconnector> initColourReconnector;
  // Describe a concrete class with persistent data.

  ColourReconnector & operator=(const ColourReconnector &);
  //  Private and non-existent assignment operator.

  int    _ClReco;
  double _PReco;

};


}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of ColourReconnector.
template <>
struct BaseClassTrait<Herwig::ColourReconnector,1> {
  typedef Pythia7::HandlerBase NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ColourReconnector>: public ClassTraitsBase<Herwig::ColourReconnector> {
  static string className() { return "/Herwig++/ColourReconnector"; }
  // Return the class name.
  static string library() { return "libHwHadronization.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "ColourReconnector.icc"

#endif /* HERWIG_ColourReconnector_H */
