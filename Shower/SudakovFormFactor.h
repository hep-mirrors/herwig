// -*- C++ -*-
#ifndef HERWIG_SudakovFormFactor_H
#define HERWIG_SudakovFormFactor_H
//
// This is the declaration of the <!id>SudakovFormFactor<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This is the abstract class from which all different types of <BR>
// Sudakov form factors derive from.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun.html">SplitFun.h</a>, <BR>
// <a href="http:ShowerAlpha.html">ShowerAlpha.h</a>. <BR>
// <a href="http:SplittingGenerator.html">SplittingGenerator.h</a>, <BR>
// <a href="http:QtoQGSudakovFormFactor.html">QtoQGSudakovFormFactor.h</a>, <BR>
// <a href="http:QtoQGammaSudakovFormFactor.html">QtoQGammaSudakovFormFactor.h</a>, <BR>
// <a href="http:GtoGGSudakovFormFactor.html">GtoGGSudakovFormFactor.h</a>, <BR>
// <a href="http:GtoQQbarSudakovFormFactor.html">GtoQQbarSudakovFormFactor.h</a>.
// 

#include "Pythia7/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Config/GlobalParameters.h"
#include "SplitFun.h"
#include "ShowerAlpha.h"


namespace Herwig {

using namespace Pythia7;


class SudakovFormFactor: public Pythia7::HandlerBase {

public:

  inline SudakovFormFactor();
  inline SudakovFormFactor(const SudakovFormFactor &);
  virtual ~SudakovFormFactor();
  // Standard ctors and dtor.

  inline SudakovFormFactor( const SplitFunPtr inputSplitFun, 
			    const tShowerAlphaPtr inputShowerAlpha,
                            const Energy inputMinScale, const Energy inputMaxScale );

  virtual Energy generateNextBranching( tPartCollHdlPtr ch, 
  					const Energy startingScale,
  					const bool reverseAngularOrder = false ) const = 0;
  // Pure virtual method, to be defined in concrete derived classes.
  // It returns the scale of the next branching; if there is no 
  // branching then it returns Energy().
  // The <!id>ch<!!id> argument is used only for Initial State branching,
  // to get access to the PDFs; the <!id>reverseOrdering<!!id> is used 
  // (when it is not equal to the default, false, value) only for 
  // Final State branching of a decaying on-shell particle. 

  virtual void setupLookupTables();
  // This virtual method is defined as empty, and it should be
  // overriden only for those derived Sudakov form factor classes
  // that use lookup tables for numerical evaluations, rather
  // than using the Monte Carlo rejection (veto) method.
  // This method is called once, during initialization, by
  // the <!class>SplittingGenerator<!!class>. 
  // General methods, usable for any type of Sudakov form factor
  // that override this method, should be provided in this
  // class in the protected session.

  inline tSplitFunPtr splitFun() const;
  // It returns the pointer to the <!class>SplitFun<!!class> object.

  inline tShowerAlphaPtr alpha() const;
  // It returns the pointer to the <!class>ShowerAlpha<!!class> object.

public:

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  //***LOOKHERE*** DEFINE HERE SOME METHODS WHICH ARE USEFUL
  //               TO COMPUTE NUMERICALLY THE SUDAKOV FORM FACTOR,
  //               IN A VERY GENERAL WAY, INDEPENDENTLY FROM THE
  //               INITIAL STATE / FINAL STATE AND FROM THE
  //               KIND AND MULTIPLICITY OF THE VERTEX.
  //               I AM NOT SURE IF THIS IS POSSIBLE: MAYBE FOR
  //               ALL 1-&GT;2 SPLITFUN, BUT NOT PROBABLY FOR 1-&GT;3.
  //               EVENTUALLY ASSUME 1-&GT;2 FOR THIS DEFAULT CASE.
  //               NOTICE THAT IS IMPORTANT TO DEFINE THESE METHODS
  //               IN THE PROTECTED PART BECAUSE THEY SHOULD BE
  //               USED ONLY BY THE DERIVED CLASSES.

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

  static AbstractClassDescription<SudakovFormFactor> initSudakovFormFactor;
  // Describe an abstract base class with persistent data.

  SudakovFormFactor & operator=(const SudakovFormFactor &);
  //  Private and non-existent assignment operator.

  SplitFunPtr _splitFun;
  tShowerAlphaPtr _alpha;
  Energy _minScale;
  Energy _maxScale;

  //***LOOKHERE*** DEFINE HERE SOME DATA STRUCTURE WHICH IS USEFUL
  //               TO IMPLEMENT THE GENERAL METHODS FOR THE NUMERIC
  //               EVALUATION WHICH ARE DEFINED IN THE PROTECTED PART

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of SudakovFormFactor.
template <>
struct BaseClassTrait<Herwig::SudakovFormFactor,1> {
  typedef Pythia7::HandlerBase NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::SudakovFormFactor>: public ClassTraitsBase<Herwig::SudakovFormFactor> {
  static string className() { return "/Herwig++/SudakovFormFactor"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "SudakovFormFactor.icc"

#endif /* HERWIG_SudakovFormFactor_H */
