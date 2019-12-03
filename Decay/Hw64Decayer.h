// -*- C++ -*-
//
// Hw64Decayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Hw64Decayer_H
#define HERWIG_Hw64Decayer_H
//
// This is the declaration of the Hw64Decayer class.
//
#include <ThePEG/Config/ThePEG.h>
#include <HwDecayerBase.h>
#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Repository/Strategy.fh>
#include <fstream>

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 * <code>Hw64Decayer</code> is a class that defines all the general routines 
 * used in HERWIG++ to imitate the HERWIG 6.4 decays. The goal is to have an exact
 * copy of HERWIG 6.4 decay routines. This will allow for easy 'callibration'
 * of the new C++ code with the old Fortran code.
 *
 *  This class handles the non-partonic decays. In general it is used for
 *  exclusive meson and baryon decays. Three different matrix elements are supported
 *
 *  - MECode=0 flat-phase space.
 *  - MECode=100 free V-A matrix element
 *  - MECode=101 bound V-A matrix element
 *
 * @see HeavyDecayer
 * @see QuarkoniumDecayer
 * @see Decayer
 * 
 */
class Hw64Decayer: public HwDecayerBase {

public:

  /**
   * Default constructor
   */
  Hw64Decayer() : MECode(0),_masstry(50) {} 

  /**
   * return true if this decayer can perfom the decay specified by the
   * given decay mode.
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;

  /**
   * for a given decay mode and a given particle instance, perform the
   * decay and return the decay products.
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

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * Weighting of phase space for V-A matrix elements
   */
  static double VAWt(Energy2 t0, Energy2 t1, Energy2 t2, InvEnergy4 t3);

  /**
   * Take an array of momenta and set the momentum member of the particles.
   * @param moms The input momenta to be assigned to the particles.
   * @param particles The particles whose momenta is to be set.
   * @param out The particles outputted with their momenta set.
   */
  void setParticleMomentum(ParticleVector & out, const cPDVector & particles, 
			   const vector<Lorentz5Momentum> & moms) const {
    unsigned int numProds = particles.size();
    for(unsigned int ix=0;ix<numProds;++ix)
      out.push_back(particles[ix]->produceParticle(moms[ix]));
  }

private:

  /**
   *  Private and non-existent assignment operator.
   */
  const Hw64Decayer & operator=(const Hw64Decayer &) = delete;

private:

  /**
   *  The code for the matrix element being used.
   */
  int MECode;

  /**
   *  Maximum number of attempts to generate the off-shell masses
   */
  unsigned int _masstry;
};

}

#endif /* HERWIG_Hw64Decayer_H */
