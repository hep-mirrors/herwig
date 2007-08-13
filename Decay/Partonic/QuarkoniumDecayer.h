// -*- C++ -*-
#ifndef HERWIG_QuarkoniumDecayer_H
#define HERWIG_QuarkoniumDecayer_H
//
// This is the declaration of the QuarkoniumDecayer class.
//

#include <ThePEG/Config/ThePEG.h>
#include <PartonicDecayerBase.h>
#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Repository/Strategy.fh>
#include <fstream>

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The QuarkoniumDecayer class is designed for the partonic decay of bottom and charmonium
 *  resonances. In general it is used for decays of the type:
 *    - \f$q,\bar{q}\f$ decay to a quark-antiquark pair generally using phase space,
 *      \e i.e. MECode=0.
 *
 *    - \f$g,g\f$ decay to two gluons normally using phase space,
 *      \e i.e. MECode=0.
 *
 *    - \f$g,g,g\f$ decay to three gluons, this will normally use the Ore-Powell 
 *      matrix element, \e i.e. MECode=130.
 *
 *    - \f$g,g,\gamma\f$ decay to two gluons and a photon, this will normally use
 *      the Ore-Powell matrix element, \e i.e. MECode=130.
 * 
 *
 *  This class supports two values of the MECode variable which can be set using
 *  the interface
 *
 *  - MECode=0   flat-phase space
 *  - MECode=130 The Ore-Powell onium matrix element.
 *
 *  This is designed to be the same as the FORTRAN HERWIG routine.
 *
 * @see HeavyDecayer
 * @see Hw64Decayer
 * @see Decayer
 *
 */
class QuarkoniumDecayer: public PartonicDecayerBase {

public:

  /**
   * Standard ctors and dtor
   */
  inline QuarkoniumDecayer();

  /**
   * Check if this decayer can perfom the decay for a particular mode
   * @param parent The decaying particle
   * @param children The decay products
   * @return true If this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(tcPDPtr parent, const PDVector & children) const;
  
  /**
   *  Perform the decay of the particle to the specified decay products
   * @param parent The decaying particle
   * @param children The decay products
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const PDVector & children) const;

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

  /**
   * Standard Persistent stream methods
   */
  void persistentOutput(PersistentOStream &) const;

  /**
   * Standard Persistent stream methods
   */
  void persistentInput(PersistentIStream &, int);

   /**
    * Standard clone methods
    */
protected:

   /**
    * Standard clone methods
    */
   inline virtual IBPtr clone() const;

   /**
    * Standard clone methods
    */
   inline virtual IBPtr fullclone() const;

private:

  /**
   * Describe a concrete class with persistant decay
   */
  static ClassDescription<QuarkoniumDecayer> initQuarkoniumDecayer;

  /**
   *  Private and non-existent assignment operator.
   */
  const QuarkoniumDecayer & operator=(const QuarkoniumDecayer &);

private:

  /**
   *  The code for the type of matrix element being used.
   */
  int MECode;
};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * QuarkoniumDecayer.
 */
template <>
struct BaseClassTrait<Herwig::QuarkoniumDecayer,1> {
  /** Typedef of the base class of QuarkoniumDecayer. */
  typedef Herwig::PartonicDecayerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * QuarkoniumDecayer class.
 */
template <>
struct ClassTraits<Herwig::QuarkoniumDecayer>: 
    public ClassTraitsBase<Herwig::QuarkoniumDecayer> {
  /** Return the class name. */
  static string className() { return "Herwig::QuarkoniumDecayer"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwPartonicDecay.so"; }
};

/** @endcond */

}

#include "QuarkoniumDecayer.icc"

#endif /* HERWIG_QuarkoniumDecayer_H */
