// -*- C++ -*-
#ifndef HERWIG_QtoQGammaSudakovFormFactor_H
#define HERWIG_QtoQGammaSudakovFormFactor_H
//
// This is the declaration of the <!id>QtoQGammaSudakovFormFactor<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This (concrete) class describes the properties of the <BR>
// Sudakov form factor for the <I>Q-&GT;QGamma</I> splitting.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SudakovFormFactor.html">SudakovFormFactor.h</a>.
// 

#include "SudakovFormFactor.h"


namespace Herwig {

using namespace Pythia7;

class QtoQGammaSudakovFormFactor: public SudakovFormFactor {

public:

  inline QtoQGammaSudakovFormFactor();
  inline QtoQGammaSudakovFormFactor(const QtoQGammaSudakovFormFactor &);
  virtual ~QtoQGammaSudakovFormFactor();
  // Standard ctors and dtor.

  inline QtoQGammaSudakovFormFactor( const SplitFunPtr inputSplitFun, 
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

  QtoQGammaSudakovFormFactor & operator=(const QtoQGammaSudakovFormFactor &);
  //  Private and non-existent assignment operator.

};

}

#include "QtoQGammaSudakovFormFactor.icc"

#endif /* HERWIG_QtoQGammaSudakovFormFactor_H */
