// -*- C++ -*-
//
// ShowerConfig.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerConfig_H
#define HERWIG_ShowerConfig_H
//
// This is the declaration of the ShowerConfig class.

#include "ThePEG/Config/ThePEG.h"
#include "Base/ShowerParticle.fh"
#include "Base/ShowerKinematics.fh" 
#include "Base/SudakovFormFactor.fh"

namespace Herwig { 
using namespace ThePEG;

/** \ingroup Shower
 *  
 *  Handy header file to be included in all Shower classes.
 *  It contains only some useful typedefs.
 */
  
  /**
   * Pointer to a ColourLine
   */
  typedef Ptr<ThePEG::ColourLine>::pointer ColinePtr;
  
  /**
   * Transient Pointer to a ColourLine
   */
  typedef Ptr<ThePEG::ColourLine>::transient_pointer tColinePtr;
  
  /**
   * A pair of ColourLine pointers
   */
  typedef pair<ColinePtr,ColinePtr> ColinePair;

  /**
   * A pair of transient ColourLine pointers
   */
  typedef pair<tColinePtr,tColinePtr> tColinePair;

  /**
   *  A Vector of ShowerParticle pointers
   */
  typedef vector<ShowerParticlePtr> ShowerParticleVector;

  /**
   *  A Vector of transient ShowerParticle pointers
   */
  typedef vector<tShowerParticlePtr> tShowerParticleVector;

  /**
   *  Definition of the IdList for branchings
   */
  typedef vector<long> IdList;

  namespace ShowerInteraction {
    /**
     *  Enum for the type of interaction
     */
    enum Type { UNDEFINED=-1, QCD, QED };
  }

  /**
   *  typedef to pair the SudakovFormFactor and the particles in a branching
   */
  typedef pair<SudakovPtr,IdList> BranchingElement;

  /**
   *  typedef to pair the PDG code of the particle and the BranchingElement
   */
  typedef multimap<long,BranchingElement> BranchingList;

  /**
   *  typedef to create a structure which can be inserted into a BranchingList
   */
  typedef pair<long, BranchingElement> BranchingInsert; 

  /** \ingroup Shower
   *  The branching struct is used to store information on the branching.
   *  The kinematics variable is a pointer to the ShowerKinematics for the branching
   *  The sudakov variable is a pointer to the SudakovFormFactor for the branching
   *  The ids  variable is the list of particles in the branching
   */
  struct Branching {
    
    /**
     *  Pointer to the ShowerKinematics object for the branching
     */
    ShoKinPtr kinematics;
    
    /**
     *  PDG codes of the particles in the branching
     */
    IdList ids; 
    
    /**
     *  The SudakovFormFactor for the branching
     */
    tSudakovPtr sudakov;
    
    /**
     *  Constructor for the struct
     * @param a pointer to the ShowerKinematics object for the branching
     * @param c PDG codes of the particles in the branching
     * @param d The SudakovFormFactor for the branching
     */
    Branching(ShoKinPtr a, IdList c,tSudakovPtr d) : kinematics(a), ids(c), sudakov(d) {}
    
    /**
     *  Default constructor
     */
    Branching() {}
  };
}

#endif // HERWIG_ShowerConfig_H 

