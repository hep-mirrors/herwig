// -*- C++ -*-
//
// DarkQuarkoniumDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DarkQuarkoniumDecayer_H
#define HERWIG_DarkQuarkoniumDecayer_H
//
// This is the declaration of the DarkQuarkoniumDecayer class.
//

#include "PartonicDecayerBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The DarkQuarkoniumDecayer class is designed for the partonic decay of bottom and charmonium
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
class DarkQuarkoniumDecayer: public PartonicDecayerBase {

public:

  /**
   * Standard ctors and dtor
   */
  DarkQuarkoniumDecayer();

  /**
   * Check if this decayer can perfom the decay for a particular mode
   * @param parent The decaying particle
   * @param children The decay products
   * @return true If this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;
  
  /**
   *  Perform the decay of the particle to the specified decay products
   * @param parent The decaying particle
   * @param children The decay products
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const tPDVector & children) const;


  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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
   virtual IBPtr clone() const;

   /**
    * Standard clone methods
    */
   virtual IBPtr fullclone() const;

private:

  /**
   *  Private and non-existent assignment operator.
   */
  const DarkQuarkoniumDecayer & operator=(const DarkQuarkoniumDecayer &) = delete;

private:

  /**
   *  The code for the type of matrix element being used.
   */
  int MECode;
};

}

#endif /* HERWIG_DarkQuarkoniumDecayer_H */
