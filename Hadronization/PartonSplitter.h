// -*- C++ -*-
#ifndef HERWIG_PartonSplitter_H
#define HERWIG_PartonSplitter_H
//
// This is the declaration of the <!id>PartonSplitter<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class does all of the nonperturbative parton splittings needed <BR>
// immediately after the end of the showering (both initial and final), <BR>
// as very first step of the cluster hadronization.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:GlobalParameters.html">GlobalParameters.h</a>.
// 

#include "CluHadConfig.h"
#include "Pythia7/Handlers/HandlerBase.h"
#include "Herwig++/Config/GlobalParameters.h"


namespace Herwig {


using namespace Pythia7;


class PartonSplitter: public Pythia7::HandlerBase {

public:

  inline PartonSplitter();
  inline PartonSplitter(const PartonSplitter &);
  virtual ~PartonSplitter();
  // Standard ctors and dtor.

public:

  void split(const tPVector & tagged, tStepPtr pstep);
  // It does the nonperturbative splitting of:
  // time-like gluons; space-like gluons; space-like sea-quark (or antiquark)
  // (in the latter case, the produced soft gluon is then splitted in valence
  // quarks). The new particles produced are added to the event record.
 
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

  static ClassDescription<PartonSplitter> initPartonSplitter;
  // Describe a concrete class with persistent data.

  PartonSplitter & operator=(const PartonSplitter &);
  //  Private and non-existent assignment operator.

  void splitTimeLikeGluon(tcPPtr ptrGluon,                   // input        
			  PPtr & ptrQ, PPtr & ptrQbar);      // output
  // Given in input a pointer to a time-like gluon, it forces
  // a nonperturbative quark - anti-quark splitting, returning
  // the pointers to these produced two new particles. 
  // If something wrong happens, it will returns null pointers.
  
  void splitSpaceLikeGluon(tcPPtr ptrGluon,                  // input       
			   PPtr & ptrQ, PPtr & ptrQbar);     // output
  // Given in input a pointer to a space-like gluon, it forces
  // a nonperturbative quark - anti-quark splitting, returning
  // the pointers to these produced two new particles. 
  // If something wrong happens, it will returns null pointers.
  
  void splitSpaceLikeSeaQuark(tcPPtr ptrSeaQ0,                   // input
  			      PPtr & ptrGluon, PPtr & ptrSeaQ1); // output
  // Given in input a pointer to a space-like sea quark, it forces 
  // a nonperturbative soft gluon emission, returning the pointers 
  // to the emitted gluon and the sea quark after the emission. 
  // If something wrong happens, it will return null pointers.

  void debuggingInfo(const tPVector & tagged, const set<tPPtr> & newPartons);
  // Print full information for debugging.

  GlobParamPtr _globalParameters;  

};


}

// CLASSDOC OFF

namespace Pythia7 {


// The following template specialization informs Pythia7 about the
// base class of PartonSplitter.
template <>
struct BaseClassTrait<Herwig::PartonSplitter,1> {
  typedef Pythia7::HandlerBase NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::PartonSplitter>:
    public ClassTraitsBase<Herwig::PartonSplitter> {
  static string className() { return "/Herwig++/PartonSplitter"; }
  // Return the class name.
  static string library() { return "libHwHadronization.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};


}

#include "PartonSplitter.icc"

#endif /* HERWIG_PartonSplitter_H */
