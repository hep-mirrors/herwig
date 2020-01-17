// -*- C++ -*-
//
// SOPHTY.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SOPHTY_H
#define HERWIG_SOPHTY_H
//
// This is the declaration of the SOPHTY class.
//

#include "DecayRadiationGenerator.h"
#include "FFDipole.fh"
#include "IFDipole.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 * Here is the documentation of the SOPHTY class.
 *
 * @see \ref SOPHTYInterfaces "The interfaces"
 * defined for SOPHTY.
 */
class SOPHTY: public DecayRadiationGenerator {

public:

  /**
   *  Default constructor
   */
  SOPHTY() : colouredOption_(0) {}

  /**
   *  Member to generate the photons in the decay. This must be implemented
   *  in classes inheriting from this one to produce the radiation.
   * @param p The decaying particle
   * @param children The decay products
   * @param decayer The decayer for with decay mode
   * @return The decay products with additional radiation
   */
  virtual ParticleVector generatePhotons(const Particle & p,ParticleVector children,
					 tDecayIntegratorPtr decayer);

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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SOPHTY & operator=(const SOPHTY &) = delete;

private:

  /**
   *  The final-final dipole
   */
  FFDipolePtr FFDipole_;

  /**
   *  The initial-final dipole
   */
  IFDipolePtr IFDipole_;

  /**
   *  Option for the treatment of radiation from coloured particles
   */
  unsigned int colouredOption_;
};

}

#endif /* HERWIG_SOPHTY_H */
