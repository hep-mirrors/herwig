#ifndef HERWIG_SextetParticles_H
#define HERWIG_SextetParticles_H

#include <ThePEG/PDT/EnumParticles.h>

namespace ThePEG {

/**
 * The ParticleID namespace defines the ParticleCodes enumeration.
 */
namespace ParticleID {

/** \ingroup Utilities
 * Enumeration to give identifiers to PDG codes for sextet particles
 */  
enum SextetCodes {
  ScalarDQSingletY43 = 6000221, 
    /**< \f$\Phi_{6,1,4/3}\f$*/
  ScalarDQSingletY13 = 6000211,
    /**< \f$\Phi_{6,1,1/3}\f$*/
  ScalarDQSingletY23 = 6000111,
    /**< \f$\Phi_{6,1,-2/3}\f$*/
  ScalarDQSingletY43bar = -6000221, 
    /**< \f$\bar\Phi_{6,1,4/3}\f$*/
  ScalarDQSingletY13bar = -6000211,
    /**< \f$\bar\Phi_{6,1,1/3}\f$*/
  ScalarDQSingletY23bar = -6000111,
    /**< \f$\bar\Phi_{6,1,-2/3}\f$*/
  ScalarDQTripletP = 6001221, 
    /**< \f$\Phi^+_{6,3,1/3}\f$*/
  ScalarDQTriplet0 = 6001211,
    /**< \f$\Phi^0_{6,3,1/3}\f$*/
  ScalarDQTripletM = 6001111,
    /**< \f$\Phi^-_{6,3,1/3}\f$*/
  ScalarDQTripletPbar = -6001221, 
    /**< \f$\bar\Phi^+_{6,3,1/3}\f$*/
  ScalarDQTriplet0bar = -6001211,
    /**< \f$\bar\Phi^0_{6,3,1/3}\f$*/
  ScalarDQTripletMbar = -6001111,
    /**< \f$\bar\Phi^-_{6,3,1/3}\f$*/
  VectorDQY16P = 6000123,
    /**< \f$V^{+\frac12\mu}_{6,2,-1/6}\f$*/
  VectorDQY16M = 6000113,
    /**< \f$V^{-\frac12\mu}_{6,2,-1/6}\f$*/
  VectorDQY16Pbar = -6000123,
    /**< \f$\bar V^{+\frac12\mu}_{6,2,-1/6}\f$*/
  VectorDQY16Mbar = -6000113,
    /**< \f$\bar V^{-\frac12\mu}_{6,2,-1/6}\f$*/
  VectorDQY56P = 6000223,
    /**< \f$V^{+\frac12\mu}_{6,2,5/6}\f$*/
  VectorDQY56M = 6000213,
    /**< \f$V^{-\frac12\mu}_{6,2,5/6}\f$*/
  VectorDQY56Pbar = -6000223,
    /**< \f$\bar V^{+\frac12\mu}_{6,2,5/6}\f$*/
  VectorDQY56Mbar = -6000213
    /**< \f$\bar V^{-\frac12\mu}_{6,2,5/6}\f$*/
};

}

}

#endif


