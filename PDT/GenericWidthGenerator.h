// -*- C++ -*-
//
// GenericWidthGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GenericWidthGenerator_H
#define HERWIG_GenericWidthGenerator_H
//
// This is the declaration of the GenericWidthGenerator class.
//
#include "ThePEG/PDT/WidthGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/DecayMode.h"
#include "GenericWidthGenerator.fh"
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Utilities/Interpolator.h"
#include "GenericMassGenerator.h"
#include <iostream>

namespace Herwig {
using namespace ThePEG;

  /**
   *  Declare ModelGenerator class as must be friend to set the particle
   */
  class ModelGenerator;

  /**
   * Typedef to define a DecayMoap
   */
  typedef Selector<tDMPtr> DecayMap;


/** \ingroup PDT
 *
 * The GenericWidthGenerator class is designed to automatically
 * calculate the running width for a given particle using information from
 * the decayModes and the Decayers to construct the running width.
 *
 * It also gives us the option of selecting the decay modes for a particle
 * based on the mass.
 *
 * @see WidthGenerator
 * @see DecayIntegrator
 * @see GenericMassGenerator
 */
class GenericWidthGenerator: public WidthGenerator {

  /**
   *  ModelGenerator class as must be friend to set the particle
   */
  friend class ModelGenerator;

public:

  /**
   * A friend class so the off-shell matrix elements can be integrated.
   */
  friend class TwoBodyAllOnCalculator;

public:

  /**
   * Default constructor
   */
  GenericWidthGenerator()
    : mass_(), prefactor_(1.), initialize_(false),output_(false),
      BRnorm_(true),npoints_(50),
      BRminimum_(0.01), intOrder_(1), twoBodyOnly_(false)
  {}

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

public:

  /**
   * Return true if this width generator can handle the given particle type.
   * @param part The particle data pointer of the particle.
   * @return True if this class can handle the particle and false otherwise
   */
  virtual bool accept(const ParticleData & part) const {
    if(!particle_) return false;
    return part.id() == particle_->id() ||
      ( part.CC() && part.CC()->id() == particle_->id() );
  }

  /** @name Members to calculate the width and decay modes. */
  //@{
  /**
   * Calculate the width.
   * @param part The particle data pointer of the particle.
   * @param m The scale for the width calculation
   * @return The width at the mass given.
   */
  virtual Energy width(const ParticleData & part, Energy m) const;

  /**
   * Initialize the given decay map for the given particle type.
   * @param part The particle data pointer of the particle.
   * @return The decay map
   */
  virtual DecayMap rate(const ParticleData & part) const {
    return part.decaySelector();
  }

  /**
   * Return a decay map for a given particle instance. This allows us to
   * vary the branching ratios as a function of the particles mass.
   * @param  part The particle instance
   * @return The decay map
   */
  virtual DecayMap rate (const Particle & part);

  /**
   * The partial width for a given mode
   * @param m The mass, or scale, for the calculation
   * @param iloc The location of the mode in the list
   * @return The partial width for the mode.
   */
  virtual Energy partialWidth(int iloc,Energy m) const;

  /**
   *  Return the total width and the sum of the partial widths for
   *  modes which are used
   */
  virtual pair<Energy,Energy> width(Energy, const ParticleData &) const;
  //@}

  /**
   * Output the initialisation info for the database
   * @param output The stream to output the information to
   * @param header output the header.
   **/
  virtual void dataBaseOutput(ofstream & output, bool header=true);

  /**
   * Given a particle type and a mass and a width of an instance of
   * that particle type, generate a life time.
   */
  virtual Length lifeTime(const ParticleData &, Energy m, Energy w) const;

protected:

  /**
   * The \f$1\to2\f$ width for on-shell particles
   * @param q The mass, or scale, for the calculation
   * @param iloc The location of the mode in the list.
   * @return The partial width.
   */
  Energy partial2BodyWidth(int iloc,Energy q) const {
    return partial2BodyWidth(iloc,q,MEmass1_[iloc],MEmass2_[iloc]);
  }

  /**
   * The \f$1\to2\f$ width for outgoing particles which can be off-shell.
   * @param iloc The location of the mode in the list.
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the first outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @return The partial width.
   */
  virtual Energy partial2BodyWidth(int iloc,Energy m0,Energy m1,Energy m2) const;

  /**
   * Perform the set up for a mode in classes inheriting from this one
   * @param mode The decay mode
   * @param decayer The decayer for the mode.
   * @param imode The number of the mode.
   */
  virtual void setupMode(tcDMPtr mode, tDecayIntegratorPtr decayer, unsigned int imode);

  /**
   *  Access to the particle dat for inheriting classes
   */
  tPDPtr particle() const {return particle_;}

  /**
   * Set the particle
   */
  void particle(tPDPtr in) {particle_ = in;}

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans)
   ;

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

  /**
   *  Matrix element code for a given mode 
   * @param imode The mode.
   */
  int MEcode(int imode) const {return MEcode_[imode];}

  /**
   *  Coupling for a given mode
   * @param imode The mode.
   */
  double MEcoupling(int imode) const {return MEcoupling_[imode];}


  /**
   *  The on-shell mass of the particle
   */
  Energy mass() const {return mass_;}

  /**
   *  Access to the decay modes
   */

  /**
   * Initialization option for use by the inheriting classes
   */
  bool initialize() const {return initialize_;}

  /**
   * Output option for use by the inheriting classes
   */
  bool output() const {return output_;}

  /**
   *  Access to the DecayModes
   */
  vector<tDMPtr> decayModes() const {return decayModes_;}
  
private:

  /**
   * Helper function for the interface
   */
  void setParticle(string);

  /**
   * Helper function for the interface
   */
  string getParticle() const;

private:

  /**
   * Private and non-existent assignment operator.
   */
  GenericWidthGenerator & operator=(const GenericWidthGenerator &) = delete;

private:

  /**
   * The pointer to the ParticleData object for the particle for this width generator.
   */
  tPDPtr particle_;

  /**
   * The decaymodes
   */
  vector<tDMPtr> decayModes_;

  /**
   *  The tags for the DecayMode s
   */
  vector<string> decayTags_;
  
  /**
   *  The minimum mass of the decaying particle for which this decay mode is possible
   */
  vector<Energy> minMass_;

  /**
   * The on-shell mass of the particle
   */
  Energy mass_;

  /**
   * Prefactor to get the on-shell width
   */
  double prefactor_;

  /**
   * The type of ME, whether it is fixed, calculated by this class or interpolation
   */
  vector<int> MEtype_;

  /**
   * The code for the matrix element
   */
  vector<int> MEcode_;

  /**
   *  Mass of the first outgoing particle for the simple \f$1\to2\f$ ME's
   */
  vector<Energy> MEmass1_;
  /**
   *  Mass of the second outgoing particle for the simple \f$1\to2\f$ ME's
   */
  vector<Energy> MEmass2_;

  /**
   * the coupling for a given me
   */
  vector<double> MEcoupling_; 

  /**
   * is this mode used for the running width
   */
  vector<bool> modeOn_;

  /**
   * storage of the massesto set up the interpolation tables
   */
  vector<Energy> interMasses_;

  /**
   * storage of the widths to set up the interpolation tables
   */
  vector<Energy> interWidths_;

  /**
   * the number of entries in the decay table for a particular mode
   */
  vector<int> noOfEntries_;

  /**
   * initialize the generator using the particle data object
   */
  bool initialize_;

  /**
   *  Output the parameters
   */
  bool output_;

  /**
   * normalise the terms so that the partial widths for an on-shell particle are correct
   */
  bool BRnorm_;

  /**
   * number of points to use for interpolation tables
   */
  int npoints_;

  /**
   * intepolators for the running width
   */
  vector<Interpolator<Energy,Energy>::Ptr> interpolators_;

  /**
   * minimum branching ratio for the inclusion in the total running width
   */
  double BRminimum_;

  /**
   *  Order of the interpolation for the tables
   */
  unsigned int intOrder_;

  /**
   *  Whether or not to only include 2 body modes in the running
   * width calculation, higher modes flat
   */
  bool twoBodyOnly_;
};

}

#endif /* HERWIG_GenericWidthGenerator_H */
