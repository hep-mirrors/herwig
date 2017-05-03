// -*- C++ -*-
//
// HwDecayerBase.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HwDecayerBase_H
#define HERWIG_HwDecayerBase_H
//
// This is the declaration of the HwDecayerBase class.
//

#include "ThePEG/PDT/Decayer.h"
#include "Herwig/Shower/Core/Base/ShowerParticle.fh"
#include "Herwig/Shower/Core/Base/ShowerProgenitor.fh"
#include "Herwig/Shower/RealEmissionProcess.fh"
#include "HwDecayerBase.fh"

namespace Herwig {

struct Branching;

using namespace ThePEG;

/**
 * The HwDecayerBase class is the base class for Decayers in Herwig. It inherits
 * from the Decayer class of ThePEG and implements additional functionality for the 
 * output of the results to the particle database and initialization of the datbase.
 *
 * It also provide the option of specifying a class based on the DecayRadiationGenerator
 * which should be used to generate QED radiation in the decay
 *
 * @see \ref HwDecayerBaseInterfaces "The interfaces"
 * defined for HwDecayerBase.

 */
class HwDecayerBase: public Decayer {

public:

  /**
   * The default constructor.
   */
  HwDecayerBase() : _initialize(false), _dbOutput(false) {}

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;
  //@}

public:

  /**
   *  Virtual members to be overridden by inheriting classes
   *  which implement hard corrections 
   */
  //@{
  
  /**
   * Type of POWHEG correction
   */
  enum POWHEGType {No, ISR, FSR, Both};


  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return No;}


  /**
   *  Has an old fashioned ME correction
   */
  virtual bool hasMECorrection() {return false;}

  /**
   *  Initialize the ME correction
   */
  virtual void initializeMECorrection(RealEmissionProcessPtr , double & ,
				      double & );

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual RealEmissionProcessPtr applyHardMatrixElementCorrection(RealEmissionProcessPtr);

  /**
   * Apply the soft matrix element correction
   * @param initial The particle from the hard process which started the 
   * shower
   * @param parent The initial particle in the current branching
   * @param br The branching struct
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr initial,
				     ShowerParticlePtr parent,
				     Branching br);

  /**
   *  Apply the POWHEG style correction
   */
  virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr);
  //@}

protected:

  /** @name Virtual functions to replaced those from the Decayer class. 
   *  This is so that the decay and accept members of this class can handle all
   *  the more complicated features of the DecayMode class
   */
  //@{
  /**
   * Check if this decayer can perfom the decay for a particular mode
   * @param parent The decaying particle
   * @param children The decay products
   * @return true If this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const = 0;
  
  /**
   *  Perform the decay of the particle to the specified decay products
   * @param parent The decaying particle
   * @param children The decay products
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const tPDVector & children) const = 0;
  //@}

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

public:

  /**
   *  Functions for the Herwig decayer
   */
  //@{
  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const = 0;

  /**
   *  Access to the initialize variable
   */
  bool initialize() const {return _initialize;}

  /**
   *  Access the database output variable
   */
  bool databaseOutput() const {return _dbOutput;}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<HwDecayerBase> initHwDecayerBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HwDecayerBase & operator=(const HwDecayerBase &);

private:

  /**
   * perform initialisation
   */
  bool _initialize;

  /**
   * Print out database  
   */
  bool _dbOutput;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HwDecayerBase. */
template <>
struct BaseClassTrait<Herwig::HwDecayerBase,1> {
  /** Typedef of the first base class of HwDecayerBase. */
  typedef Decayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HwDecayerBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HwDecayerBase>
  : public ClassTraitsBase<Herwig::HwDecayerBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HwDecayerBase"; }
};

/** @endcond */

}

#endif /* HERWIG_HwDecayerBase_H */
