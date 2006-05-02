// -*- C++ -*-
#ifndef HERWIG_ShowerConfig_H
#define HERWIG_ShowerConfig_H
//
// This is the declaration of the ShowerConfig class.

#include "Herwig++/Config/Herwig.h"
#include "Kinematics/ShowerParticle.fh"
#include "ShowerVariables.fh"
#include "Kinematics/ShowerKinematics.fh" 

namespace Herwig { 
using namespace ThePEG;

  /** \ingroup Shower
   *  
   *  Handy header file to be included in all Shower classes.
   *  It contains only some useful typedefs.
   */
  class ThePEG::ColourLine;

  typedef Ptr<ThePEG::ColourLine>::pointer ColinePtr;
  typedef Ptr<ThePEG::ColourLine>::transient_pointer tColinePtr;
  typedef pair<ColinePtr,ColinePtr> ColinePair;
  typedef pair<tColinePtr,tColinePtr> tColinePair;

  typedef vector<ShowerParticlePtr> ShowerParticleVector;
  typedef vector<tShowerParticlePtr> tShowerParticleVector;

  typedef vector<ShoKinPtr> CollecShoKinPtr;
  typedef vector<tShoKinPtr> tCollecShoKinPtr;

  typedef vector<long> IdList;

} // end Herwig namespace

#endif // HERWIG_ShowerConfig_H 

