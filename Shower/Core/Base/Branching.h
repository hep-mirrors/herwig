// -*- C++ -*-
#ifndef HERWIG_Branching_H
#define HERWIG_Branching_H
//
// This is the declaration of the Branching struct.
//
#include "Herwig/Shower/Core/ShowerConfig.h"
#include "Herwig/Shower/Core/Base/ShowerKinematics.h"

namespace Herwig {

using namespace ThePEG;

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
     *  The type of radiation line
     */
    ShowerPartnerType type;

    /**
     *  Whether or not it keep from forced hard emisson
     */
    bool hard;

    /**
     *  Which of the children is same as incoming
     */
    unsigned int iout;
    
    /**
     *  Constructor for the struct
     * @param a pointer to the ShowerKinematics object for the branching
     * @param c PDG codes of the particles in the branching
     * @param d The SudakovFormFactor for the branching
     */
    Branching(ShoKinPtr a, IdList c,tSudakovPtr d,ShowerPartnerType t) 
      : kinematics(a), ids(c), sudakov(d), type(t), hard(false), iout(0) {}
    
    /**
     *  Default constructor
     */
    Branching() : hard(false), iout(0) {}
  };

}

#endif /* HERWIG_Branching_H */
