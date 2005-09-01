// -*- C++ -*-
#ifndef HERWIG_PartonSplitter_H
#define HERWIG_PartonSplitter_H

#include "CluHadConfig.h"
#include <ThePEG/Handlers/HandlerBase.h>
#include "Herwig++/Utilities/GlobalParameters.h"


namespace Herwig {


using namespace ThePEG;


/** \ingroup Hadronization
 *  \class PartonSplitter
 *  \brief This class splits the gluons from the end of the shower.
 *  \author Philip Stephens
 *  \author Alberto Ribon
 * 
 *  This class does all of the nonperturbative parton splittings needed 
 *  immediately after the end of the showering (both initial and final),
 *  as very first step of the cluster hadronization.
 *
 *  See also:
 *  GlobalParameters.h.
 */
class PartonSplitter: public ThePEG::HandlerBase {

public:

  /**
   * Standard ctors and dtor.
   */
  inline PartonSplitter();
  inline PartonSplitter(const PartonSplitter &);
  virtual ~PartonSplitter();

public:

  /**
   * This method does the nonperturbative splitting of:
   * time-like gluons. At the end of the shower the gluons should be
   * on a "physical" mass shell and should therefore be time-like.
   * @param tagged The tagged particles to be split
   * @param pstep Pointer to the step
   * @return The particles which were not split and the products of splitting.
   */
  tPVector split(const tPVector & tagged, tStepPtr pstep);
 
public:

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   * Standard clone methods.
   */
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;

protected:

  /**
   * Standard Interfaced virtual functions.
   */
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();

  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<PartonSplitter> initPartonSplitter;

  /**
   * Private and non-existent assignment operator.
   */
  PartonSplitter & operator=(const PartonSplitter &);

  /**
   * Given in input a pointer to a time-like gluon, it forces
   * a nonperturbative quark - anti-quark splitting, returning
   * the pointers to these produced two new particles. 
   * If something wrong happens, it will returns null pointers.
   */
  void splitTimeLikeGluon(tcPPtr ptrGluon,                   // input        
			  PPtr & ptrQ, PPtr & ptrQbar);      // output
  
  /**
   * Given in input a pointer to a space-like gluon, it forces
   * a nonperturbative quark - anti-quark splitting, returning
   * the pointers to these produced two new particles. 
   * If something wrong happens, it will returns null pointers.
   */
  //void splitSpaceLikeGluon(tcPPtr ptrGluon,                  // input       
  //         		     PPtr & ptrQ, PPtr & ptrQbar);     // output
  
  /**
   * Given in input a pointer to a space-like sea quark, it forces 
   * a nonperturbative soft gluon emission, returning the pointers 
   * to the emitted gluon and the sea quark after the emission. 
   * If something wrong happens, it will return null pointers.
   */
  //void splitSpaceLikeSeaQuark(tcPPtr ptrSeaQ0,                   // input
  //			        PPtr & ptrGluon, PPtr & ptrSeaQ1); // output

  /**
   * Print full information for debugging.
   */
  void debuggingInfo(const tPVector & tagged, const set<tPPtr> & newPartons);

  /**
   * Pointer to a Herwig::GlobalParameters object for using global variables.
   */
  GlobParamPtr _globalParameters;  

};


}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {


/**
 * The following template specialization informs ThePEG about the
 * base class of PartonSplitter.
 */
template <>
struct BaseClassTrait<Herwig::PartonSplitter,1> {
  typedef ThePEG::HandlerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::PartonSplitter>:
    public ClassTraitsBase<Herwig::PartonSplitter> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/PartonSplitter"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwHadronization.so"; }
};

}

#endif // DOXYGEN

#include "PartonSplitter.icc"

#endif /* HERWIG_PartonSplitter_H */
