// -*- C++ -*-
#ifndef HERWIG_QtoQGSudakovFormFactor_H
#define HERWIG_QtoQGSudakovFormFactor_H
//
// This is the declaration of the <!id>QtoQGSudakovFormFactor<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This (concrete) class describes the properties of the <BR>
// Sudakov form factor for the <I>Q-&GT;QG</I> splitting.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SudakovFormFactor.html">SudakovFormFactor.h</a>.
// 

#include "SudakovFormFactor.h"


namespace Herwig {

using namespace Pythia7;

class QtoQGSudakovFormFactor: public SudakovFormFactor {

public:

  inline QtoQGSudakovFormFactor();
  inline QtoQGSudakovFormFactor(const QtoQGSudakovFormFactor &);
  virtual ~QtoQGSudakovFormFactor();
  // Standard ctors and dtor.

  inline QtoQGSudakovFormFactor( const SplitFunPtr inputSplitFun, 
				 const tShowerAlphaPtr inputShowerAlpha,
				 const Energy inputInfScale, const Energy inputSupScale );

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

  QtoQGSudakovFormFactor & operator=(const QtoQGSudakovFormFactor &);
  //  Private and non-existent assignment operator.

};

}

#include "QtoQGSudakovFormFactor.icc"

#endif /* HERWIG_QtoQGSudakovFormFactor_H */
