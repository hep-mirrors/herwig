// -*- C++ -*-
#ifndef HERWIG_ShowerConfig_H
#define HERWIG_ShowerConfig_H
//
// This is the declaration of the ShowerConfig class.

#include "Herwig++/Config/Herwig.h"
#include "ThePEG/Config/Complex.h" 
#include "ShowerVariables.fh"
#include "ShowerKinematics.fh"


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

  typedef vector<ShoKinPtr> CollecShoKinPtr;
  typedef vector<tShoKinPtr> tCollecShoKinPtr;

  typedef vector<long> IdList;

} // end Herwig namespace


#endif // HERWIG_ShowerConfig_H 
