// -*- C++ -*-
#ifndef HERWIG_EvtGenRandom_H
#define HERWIG_EvtGenRandom_H
//
// This is the declaration of the EvtGenRandom class.
//
#include <EvtGenBase/EvtRandomEngine.hh>
#include <ThePEG/Repository/RandomGenerator.h>

namespace Herwig {

using namespace ThePEG;

/**
 * The EvtGenRamdom class is a wrapper around the RandomGenerator class of ThePEG
 * so that when running Herwig++ with EvtGen the same random number generator is used
 * for both.
 */
class EvtGenRandom : public EvtRandomEngine {

public:

  /**
   * The Constructor.
   * @param rand Pointer to the random number generator
   */
  inline EvtGenRandom(Ptr<RandomGenerator>::pointer rand){_rand=rand;}

  /**
   *  Member to return the random number
   */
  inline double random(){return _rand->rnd();}

  inline virtual ~EvtGenRandom(){;}

private:

  /**
   *  The pointer to the random number generator
   */
  Ptr<RandomGenerator>::pointer _rand;
};
}
#endif /* HERWIG_EvtGenRandom_H */
