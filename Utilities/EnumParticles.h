#ifndef HERWIG_EnumParticles_H
#define HERWIG_EnumParticles_H

#include <ThePEG/PDT/EnumParticles.h>

namespace Herwig {

/**
 * The ParticleID namespace defines the ParticleCodes enumeration.
 */
namespace ExtraParticleID {

/** \ingroup Utilities
 * Enumeration to give identifiers to PDG codes for special Herwig++ particles
 */  
enum ParticleCodes {
  Cluster = 81, 
    /**< Herwig++ Cluster for Hadronization*/
  Remnant = 82
    /**< Herwig++ beam Remnant */
};

}

}

#endif


