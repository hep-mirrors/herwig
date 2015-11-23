// -*- C++ -*-
#ifndef HERWIG_Branching_H
#define HERWIG_Branching_H
//
// This is the declaration of the Branching struct.
//
#include "Herwig/Shower/ShowerConfig.h"
#include "Herwig/Shower/Base/ShowerKinematics.h"

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
    ShowerPartnerType::Type type;
    
    /**
     *  Constructor for the struct
     * @param a pointer to the ShowerKinematics object for the branching
     * @param c PDG codes of the particles in the branching
     * @param d The SudakovFormFactor for the branching
     */
    Branching(ShoKinPtr a, IdList c,tSudakovPtr d,ShowerPartnerType::Type t) 
      : kinematics(a), ids(c), sudakov(d), type(t) {}
    
    /**
     *  Default constructor
     */
    Branching() {}
  };

}

#endif /* HERWIG_Branching_H */
