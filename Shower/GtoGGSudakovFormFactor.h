// -*- C++ -*-
#ifndef HERWIG_GtoGGSudakovFormFactor_H
#define HERWIG_GtoGGSudakovFormFactor_H
//
// This is the declaration of the <!id>GtoGGSudakovFormFactor<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This (concrete) class describes the properties of the <BR>
// Sudakov form factor for the <I>G-&GT;GG</I> splitting.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SudakovFormFactor.html">SudakovFormFactor.h</a>.
// 

#include "SudakovFormFactor.h"


namespace Herwig {

using namespace Pythia7;

class GtoGGSudakovFormFactor: public SudakovFormFactor {

public:

  inline GtoGGSudakovFormFactor();
  inline GtoGGSudakovFormFactor(const GtoGGSudakovFormFactor &);
  virtual ~GtoGGSudakovFormFactor();
  // Standard ctors and dtor.

  inline GtoGGSudakovFormFactor( const SplitFunPtr inputSplitFun, 
				 const tShowerAlphaPtr inputShowerAlpha,
				 const Energy inputInfScale, const Energy inputSupScale );

  virtual Energy generateNextBranching( tPartCollHdlPtr ch, 
  					const Energy startingScale,
  					const bool reverseAngularOrder = false ) const;
  // It returns the scale of the next branching; if there is no 
  // branching then it returns Energy().
  // The <!id>ch<!!id> argument is used only for Initial State branching,
  // to get access to the PDFs; the <!id>reverseOrdering<!!id> is used 
  // (when it is not equal to the default, false, value) only for 
  // Final State branching of a decaying on-shell particle. 

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

private:

  static ClassDescription<GtoGGSudakovFormFactor> initGtoGGSudakovFormFactor;
  // Describe an concrete class with persistent data.

  GtoGGSudakovFormFactor & operator=(const GtoGGSudakovFormFactor &);
  //  Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of GtoGGSudakovFormFactor.
template <>
struct BaseClassTrait<Herwig::GtoGGSudakovFormFactor,1> {
  typedef Herwig::SudakovFormFactor NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::GtoGGSudakovFormFactor>: public ClassTraitsBase<Herwig::GtoGGSudakovFormFactor> {
  static string className() { return "/Herwig++/GtoGGSudakovFormFactor"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "GtoGGSudakovFormFactor.icc"

#endif /* HERWIG_GtoGGSudakovFormFactor_H */
