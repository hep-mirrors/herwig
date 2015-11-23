#ifndef HERWIG_EnumParticles_H
#define HERWIG_EnumParticles_H

#include <ThePEG/PDT/EnumParticles.h>

namespace ThePEG {

/**
 * The ParticleID namespace defines the ParticleCodes enumeration.
 */
namespace ParticleID {

/** \ingroup Utilities
 * Enumeration to give identifiers to PDG codes for special Herwig particles
 */  
enum HerwigCodes {
  Cluster = 81, 
    /**< Herwig Cluster for Hadronization*/
  Remnant = 82
    /**< Herwig beam Remnant */
};

}

}

#endif


