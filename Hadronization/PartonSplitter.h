// -*- C++ -*-
#ifndef HERWIG_PartonSplitter_H
#define HERWIG_PartonSplitter_H
/*! \class Herwig::PartonSplitter PartonSplitter.h "Herwig++\Hadronization\PartonSplitter.h"
 * \brief This class splits the gluons from the end of the shower.
 * \author Philip Stephens
 * \author Alberto Ribon
 * \ingroup Hadronization
 *
 * This class does all of the nonperturbative parton splittings needed 
 * immediately after the end of the showering (both initial and final),
 * as very first step of the cluster hadronization.
 *
 * See also:
 * GlobalParameters.h.
 */

#include "CluHadConfig.h"
#include <ThePEG/Handlers/HandlerBase.h>
#include "Herwig++/Utilities/GlobalParameters.h"


namespace Herwig {


using namespace ThePEG;


class PartonSplitter: public ThePEG::HandlerBase {

public:

  inline PartonSplitter();
  inline PartonSplitter(const PartonSplitter &);
  virtual ~PartonSplitter();
  // Standard ctors and dtor.

public:

  void split(const tPVector & tagged, tStepPtr pstep);
  /*!< This method does the nonperturbative splitting of:
   * time-like gluons. At the end of the shower the gluons should be
   * on a "physical" mass shell and should therefore be time-like.
   */
 
public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  //!< Standard Init function used to initialize the interfaces.

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
  //!< Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  //!< Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<PartonSplitter> initPartonSplitter;
  //!< Describe a concrete class with persistent data.

  PartonSplitter & operator=(const PartonSplitter &);
  //!<  Private and non-existent assignment operator.

  void splitTimeLikeGluon(tcPPtr ptrGluon,                   // input        
			  PPtr & ptrQ, PPtr & ptrQbar);      // output
  /*!< Given in input a pointer to a time-like gluon, it forces
   * a nonperturbative quark - anti-quark splitting, returning
   * the pointers to these produced two new particles. 
   * If something wrong happens, it will returns null pointers.
   */
  
  //void splitSpaceLikeGluon(tcPPtr ptrGluon,                  // input       
  //			   PPtr & ptrQ, PPtr & ptrQbar);     // output
  // Given in input a pointer to a space-like gluon, it forces
  // a nonperturbative quark - anti-quark splitting, returning
  // the pointers to these produced two new particles. 
  // If something wrong happens, it will returns null pointers.
  
  //void splitSpaceLikeSeaQuark(tcPPtr ptrSeaQ0,                   // input
  //			      PPtr & ptrGluon, PPtr & ptrSeaQ1); // output
  // Given in input a pointer to a space-like sea quark, it forces 
  // a nonperturbative soft gluon emission, returning the pointers 
  // to the emitted gluon and the sea quark after the emission. 
  // If something wrong happens, it will return null pointers.

  void debuggingInfo(const tPVector & tagged, const set<tPPtr> & newPartons);
  //!< Print full information for debugging.

  GlobParamPtr _globalParameters;  
  //!< Pointer to a Herwig::GlobalParameters object for using global variables.

};


}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {


// The following template specialization informs ThePEG about the
// base class of PartonSplitter.
template <>
struct BaseClassTrait<Herwig::PartonSplitter,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
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

#endif // DOXYGEN

#include "PartonSplitter.icc"

#endif /* HERWIG_PartonSplitter_H */
