// -*- C++ -*-
#ifndef HERWIG_ShowerConfig_H
#define HERWIG_ShowerConfig_H
//
// This is the declaration of the ShowerConfig class.

#include "Herwig++/Config/Herwig.h"
#include "ThePEG/Config/Complex.h" 


namespace Herwig { 

  using namespace ThePEG;

  /** \ingroup Shower
   *  
   *  Handy header file to be included in all Shower classes.
   *  It contains only some useful typedefs.
   */

  class ThePEG::ColourLine;
  typedef Ptr<ThePEG::ColourLine>::pointer ShoColinePtr;
  typedef Ptr<ThePEG::ColourLine>::transient_pointer tShoColinePtr;
  typedef pair<ShoColinePtr,ShoColinePtr> ShoColinePair;
  typedef pair<tShoColinePtr,tShoColinePtr> tShoColinePair;

  class ShowerParticle;
  typedef Ptr<ShowerParticle>::pointer ShowerParticlePtr;
  typedef Ptr<ShowerParticle>::transient_pointer tShowerParticlePtr;
  typedef vector<ShowerParticlePtr> ShowerParticleVector;
  typedef vector<tShowerParticlePtr> tShowerParticleVector;

  class SplittingFunction;
  typedef Ptr<SplittingFunction>::pointer SplittingFnPtr;
  typedef Ptr<SplittingFunction>::transient_pointer tSplittingFnPtr;

  class SudakovFormFactor;
  typedef Ptr<SudakovFormFactor>::pointer SudakovPtr;
  typedef Ptr<SudakovFormFactor>::transient_pointer tSudakovPtr;

  class ShowerKinematics;
  typedef Ptr<ShowerKinematics>::pointer ShoKinPtr;
  typedef Ptr<ShowerKinematics>::transient_pointer tShoKinPtr;
  typedef vector<ShoKinPtr> CollecShoKinPtr;
  typedef vector<tShoKinPtr> tCollecShoKinPtr;

  // For rho-D (spin density matrix and decay matrix) matrices  
  // // // typedef std::complex<double> Complex; // Now in ThePEG/Config/Complex.h 
  typedef vector<Complex>       ComplexVector;
  typedef vector<ComplexVector> ComplexMatrix;

  //class MECorrection;
  //typedef Ptr<MECorrection>::pointer MECorrectionPtr;
  //typedef Ptr<MECorrection>::transient_pointer tMECorrectionPtr;

  class ShowerVariables;
  typedef Ptr<ShowerVariables>::pointer ShowerVarsPtr;
  typedef Ptr<ShowerVariables>::transient_pointer tShowerVarsPtr;

  class ShowerAlpha;
  typedef Ptr<ShowerAlpha>::pointer ShowerAlphaPtr;
  typedef Ptr<ShowerAlpha>::transient_pointer tShowerAlphaPtr;

  typedef vector<long> IdList;

} // end Herwig namespace


#endif // HERWIG_ShowerConfig_H 
