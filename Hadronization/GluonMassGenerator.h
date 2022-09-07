// -*- C++ -*-
//
// GluonMassGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_GluonMassGenerator_H
#define Herwig_GluonMassGenerator_H
//
// This is the declaration of the GluonMassGenerator class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/EventRecord/Particle.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Hadronization
 * \brief Dynamic gluon mass generator; the default returns a constant mass.
 *
 * @see \ref GluonMassGeneratorInterfaces "The interfaces"
 * defined for GluonMassGenerator.
 */
class GluonMassGenerator: public HandlerBase {

public:

  /**
   * Generate a single gluon mass with possible reference to a hard
   * scale Q and up to a maximum value
   */
  virtual Energy generate(Energy, Energy) const {
    return generate();
  }

  /**
   * Generate a single gluon mass with possible reference to a hard
   * scale Q
   */
  virtual Energy generate(Energy) const {
    return generate();
  }				

  /**
   * Generate a single gluon mass without further constraints
   */
  virtual Energy generate() const {
    return getParticleData(ThePEG::ParticleID::g)->constituentMass();
  }

  /**
   * Generate a list of n gluon masses, with a maximum available energy
   */
  list<Energy> generateMany(size_t n, Energy QMax) const {
    list<Energy> res;
    Energy m0, mu, md, ms, mg, mgmax, summg;

    mu=getParticleData(ThePEG::ParticleID::u)->constituentMass();
    md=getParticleData(ThePEG::ParticleID::d)->constituentMass();
    ms=getParticleData(ThePEG::ParticleID::s)->constituentMass();

    m0=md;
    if(mu<m0){m0=mu;}
    if(ms<m0){m0=ms;}

    if( QMax<2.0*m0*n ){
      throw Exception() << "cannot reshuffle to constituent mass shells" << Exception::eventerror;
    }

    bool repeat=true;

    while( repeat ){
      repeat=false;
      summg = 0.0*GeV;
      res.clear();
      for( size_t k = 0; k < n; ++k ){ 
        mg = generate();
        res.push_back(mg);
        summg += mg;
        if( summg > QMax - 2.0*m0*(n-k-1) ){
          repeat=true;
          break;
        }
      }
    }

    return res;
    
  }

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GluonMassGenerator & operator=(const GluonMassGenerator &) = delete;

};

}

#endif /* Herwig_GluonMassGenerator_H */
