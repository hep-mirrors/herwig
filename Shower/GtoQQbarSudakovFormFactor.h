// -*- C++ -*-
#ifndef HERWIG_GtoQQbarSudakovFormFactor_H
#define HERWIG_GtoQQbarSudakovFormFactor_H
//
// This is the declaration of the <!id>GtoQQbarSudakovFormFactor<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This (concrete) class describes the properties of the <BR>
// Sudakov form factor for the <I>G-&GT;QQbar</I> splitting.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SudakovFormFactor.html">SudakovFormFactor.h</a>.
// 

#include "SudakovFormFactor.h"


namespace Herwig {

using namespace Pythia7;

class GtoQQbarSudakovFormFactor: public SudakovFormFactor {

public:

  inline GtoQQbarSudakovFormFactor();
  inline GtoQQbarSudakovFormFactor(const GtoQQbarSudakovFormFactor &);
  virtual ~GtoQQbarSudakovFormFactor();
  // Standard ctors and dtor.

  inline GtoQQbarSudakovFormFactor( const SplitFunPtr inputSplitFun, 
				    const tShowerAlphaPtr inputShowerAlpha,
				    const Energy inputInfScale, const Energy inputSupScale, 
				    const Energy inQ0 );

  virtual Energy generateNextBranching( tPartCollHdlPtr ch, 
  					const Energy startingScale,
  					const bool reverseAngularOrder = false );
  // It returns the scale of the next branching; if there is no 
  // branching then it returns Energy().
  // The <!id>ch<!!id> argument is used only for Initial State branching,
  // to get access to the PDFs; the <!id>reverseOrdering<!!id> is used 
  // (when it is not equal to the default, false, value) only for 
  // Final State branching of a decaying on-shell particle. 

private:

  GtoQQbarSudakovFormFactor & operator=(const GtoQQbarSudakovFormFactor &);
  //  Private and non-existent assignment operator.

};

}

#include "GtoQQbarSudakovFormFactor.icc"

#endif /* HERWIG_GtoQQbarSudakovFormFactor_H */
