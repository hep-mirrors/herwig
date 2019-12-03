// -*- C++ -*-
//
// ShowerConfig.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerConfig_H
#define HERWIG_ShowerConfig_H
//
// This is the declaration of the ShowerConfig class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/PDT/ParticleData.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.fh"
#include "Herwig/Shower/QTilde/SplittingFunctions/SudakovFormFactor.fh"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Shower/ShowerInteraction.h"

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
  typedef vector<tcPDPtr> IdList;
  
  inline ShowerInteraction convertInteraction(ShowerPartnerType partner) {
    if (partner==ShowerPartnerType::QCDColourLine ||
        partner==ShowerPartnerType::QCDAntiColourLine)
      return ShowerInteraction::QCD;
    else if (partner==ShowerPartnerType::QED)
      return ShowerInteraction::QED;
    else
      return ShowerInteraction::UNDEFINED;
  }

  /**
   *  typedef to pair the SudakovFormFactor and the particles in a branching
   */
  struct BranchingElement {
    
    /**
     *  Constructor
     */ 
    BranchingElement();

    /**
     *  Constructor
     */ 
    BranchingElement(SudakovPtr sud, IdList part);

    /**
     * Destructor
     */
    ~BranchingElement();

    /**
    *   Access to the Sudakov
    */
    SudakovPtr sudakov;

    /**
     *   Access to the particles
     */
    IdList particles;

    /**
     *  Access to the charge conjugate particles
     */
    IdList conjugateParticles;
    
    /**
     * Compare two BranchingElements
     **/
    bool operator == (const BranchingElement& x) const{
      return
      sudakov == x.sudakov &&
      particles == x.particles &&
      conjugateParticles == x.conjugateParticles;
    }
  };

  /**
   *  typedef to pair the PDG code of the particle and the BranchingElement
   */
  typedef multimap<long,BranchingElement> BranchingList;

  /**
   *  typedef to create a structure which can be inserted into a BranchingList
   */
  typedef pair<long, BranchingElement> BranchingInsert;

}

namespace ThePEG {

  /** 
   * Output operator to allow the structure
   */
  PersistentOStream & operator << (PersistentOStream & os, 
				   const Herwig::BranchingElement  & x);

  /** 
   * Input operator to allow the structure
   */
  PersistentIStream & operator >> (PersistentIStream & is, 
				   Herwig::BranchingElement  & x);

}

#endif // HERWIG_ShowerConfig_H 

