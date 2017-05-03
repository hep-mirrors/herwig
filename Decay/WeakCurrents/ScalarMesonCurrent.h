// -*- C++ -*-
//
// ScalarMesonCurrent.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ScalarMesonCurrent_H
#define HERWIG_ScalarMesonCurrent_H
// This is the declaration of the ScalarMesonCurrent class.

#include "WeakDecayCurrent.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The weak current for the production of one (pseudo)-scalar meson.
 *
 *  In this case the current is given by
 *  \f[J^\mu = f_Pp_P^\mu,\f]
 *  where
 * - \f$f_P\f$ is the decay constant for the meson,
 * - \f$p_P\f$ is the momentum of the meson.
 *
 *  The outgoing mesons and their decay constants can be specified using the
 *  interfaces.
 *
 * @see WeakDecayCurrent.
 * 
 */
class ScalarMesonCurrent: public WeakDecayCurrent {

public:

  /**
   * Default constructor
   */
  ScalarMesonCurrent();

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

public:

  /** @name Methods for the construction of the phase space integrator. */
  //@{
  /**
   * Complete the construction of the decay mode for integration.
   * This version just adds the meson as the daughter of the last
   * resonance in the phase space channel.
   * @param icharge The total charge of the outgoing particles in the current.
   * @param imode   The mode in the current being asked for.
   * @param mode    The phase space mode for the integration
   * @param iloc    The location of the of the first particle from the current in
   *                the list of outgoing particles.
   * @param ires    The location of the first intermediate for the current.
   * @param phase   The prototype phase space channel for the integration.
   * @param upp     The maximum possible mass the particles in the current are
   *                allowed to have.
   * @return Whether the current was sucessfully constructed.
   */
  virtual bool createMode(int icharge,unsigned int imode,DecayPhaseSpaceModePtr mode,
			  unsigned int iloc,unsigned int ires,
			  DecayPhaseSpaceChannelPtr phase,Energy upp);

  /**
   * The particles produced by the current. This just returns the pseudoscalar
   * meson.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual tPDVector particles(int icharge, unsigned int imode, int iq, int ia);
  //@}

  /**
   * Hadronic current. This version returns the hadronic current described above.
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @param meopt Option for the calculation of the matrix element
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVectorE> 
  current(const int imode,const int ichan, Energy & scale, 
	  const ParticleVector & decay, DecayIntegrator::MEOption meopt) const;

  /**
   * Accept the decay. Checks the meson against the list
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Return the decay mode number for a given set of particles in the current. 
   * Checks the meson against the list
   * @param id The id's of the particles in the current.
   * @return The number of the mode
   */
  virtual unsigned int decayMode(vector<int> id);

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

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
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ScalarMesonCurrent> initScalarMesonCurrent;

  /**
   * Private and non-existent assignment operator.
   */
  ScalarMesonCurrent & operator=(const ScalarMesonCurrent &);

private:

  /**
   * the pdg code for the meson
   */
  vector<long> _id;

  /**
   * the decay constant
   */
  vector<Energy> _decay_constant;

  /**
   * The \f$\eta-\eta'\f$ mixing angle 
   */
  double _thetaeta;

  /**
   * The inital size of the arrays
   */
  unsigned int _initsize;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of ScalarMesonCurrent.
 */
template <>
 struct BaseClassTrait<Herwig::ScalarMesonCurrent,1> {
  /** Typedef of the base class of ScalarMesonCurrent. */
  typedef Herwig::WeakDecayCurrent NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ScalarMesonCurrent>
  : public ClassTraitsBase<Herwig::ScalarMesonCurrent> {
  /** Return the class name. */
  static string className() { return "Herwig::ScalarMesonCurrent"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwWeakCurrents.so"; }

};

/** @endcond */

}

#endif /* HERWIG_ScalarMesonCurrent_H */
