// -*- C++ -*-
#ifndef HERWIG_ShowerConfig_H
#define HERWIG_ShowerConfig_H
//
// This is the declaration of the <!id>ShowerConfig<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// Handy header file to be included in all Shower classes. <BR>
// It contains only some useful typedefs.
//

#include "Herwig++/Config/Herwig.h"
#include "Pythia7/Config/Complex.h" 


namespace Herwig { 

  using namespace Pythia7;

  class ShowerColourLine;
  typedef Ptr<ShowerColourLine>::pointer ShoColinePtr;
  typedef Ptr<ShowerColourLine>::transient_pointer tShoColinePtr;
  typedef pair<ShoColinePtr,ShoColinePtr> ShoColinePair;
  typedef pair<tShoColinePtr,tShoColinePtr> tShoColinePair;

  class ShowerParticle;
  typedef Ptr<ShowerParticle>::pointer ShoParPtr;
  typedef Ptr<ShowerParticle>::transient_pointer tShoParPtr;
  typedef vector<ShoParPtr> CollecShoParPtr;
  typedef vector<tShoParPtr> tCollecShoParPtr;

  class SplitFun;
  typedef Ptr<SplitFun>::pointer SplitFunPtr;
  typedef Ptr<SplitFun>::transient_pointer tSplitFunPtr;

  class SplitFun1to2;
  typedef Ptr<SplitFun1to2>::pointer SplitFun1to2Ptr;
  typedef Ptr<SplitFun1to2>::transient_pointer tSplitFun1to2Ptr;

  class SudakovFormFactor;
  typedef Ptr<SudakovFormFactor>::pointer SudakovFormFactorPtr;
  typedef Ptr<SudakovFormFactor>::transient_pointer tSudakovFormFactorPtr;

  class ShowerKinematics;
  typedef Ptr<ShowerKinematics>::pointer ShoKinPtr;
  typedef Ptr<ShowerKinematics>::transient_pointer tShoKinPtr;
  typedef vector<ShoKinPtr> CollecShoKinPtr;
  typedef vector<tShoKinPtr> tCollecShoKinPtr;

  // For rho-D (spin density matrix and decay matrix) matrices  
  // // // typedef std::complex<double> Complex; // Now in Pythia7/Config/Complex.h 
  typedef vector<Complex>       ComplexVector;
  typedef vector<ComplexVector> ComplexMatrix;

  class MECorrection;
  typedef Ptr<MECorrection>::pointer MECorrectionPtr;
  typedef Ptr<MECorrection>::transient_pointer tMECorrectionPtr;

  class ShowerConstrainer;
  typedef Ptr<ShowerConstrainer>::pointer ShoConstrPtr;
  typedef Ptr<ShowerConstrainer>::transient_pointer tShoConstrPtr;

  class ShowerAlpha;
  typedef Ptr<ShowerAlpha>::pointer ShowerAlphaPtr;
  typedef Ptr<ShowerAlpha>::transient_pointer tShowerAlphaPtr;

} // end Herwig namespace


#endif // HERWIG_ShowerConfig_H 
