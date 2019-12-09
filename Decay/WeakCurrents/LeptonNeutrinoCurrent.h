// -*- C++ -*-
//
// LeptonNeutrinoCurrent.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_LeptonNeutrinoCurrent_H
#define HERWIG_LeptonNeutrinoCurrent_H
//
// This is the declaration of the LeptonNeutrinoCurrent class.
//
#include "WeakCurrent.h"
#include "LeptonNeutrinoCurrent.fh"
#include "ThePEG/Helicity/LorentzSpinorBar.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;

/** \ingroup Decay
 *
 *  This class implements the weak decay current for a lepton and a neutrino.
 *  In this case the current is given by
 *  \f$J^\mu = \bar{u}(p_\nu)\gamma^\mu(1-\gamma_5)u(p_\ell)\f$ 
 *  where
 *  - \f$p_\nu\f$  is the momentum of the neutrino,
 *  - \f$p_\ell\f$ is the momentum of the charged lepton. 
 *
 * @see WeakCurrent.
 * 
 */
class LeptonNeutrinoCurrent: public WeakCurrent {

public:

  /**
   * Default constructor
   */
  LeptonNeutrinoCurrent() {
    // set up the modes in the base class
    addDecayMode(11,-12);
    addDecayMode(13,-15);
    addDecayMode(15,-16);
    setInitialModes(3);
  }

public:

  /** @name Methods for the construction of the phase space integrator. */
  //@{
  /**
   * Complete the construction of the decay mode for integration.classes inheriting
   * from this one.
   * This method is purely virtual and must be implemented in the classes inheriting
   * from WeakCurrent.
   * @param icharge   The total charge of the outgoing particles in the current.
   * @param resonance If specified only include terms with this particle
   * @param flavour Information on the required flavours of the quarks
   * @param imode     The mode in the current being asked for.
   * @param mode      The phase space mode for the integration
   * @param iloc      The location of the of the first particle from the current in
   *                  the list of outgoing particles.
   * @param ires      The location of the first intermediate for the current.
   * @param phase     The prototype phase space channel for the integration.
   * @param upp       The maximum possible mass the particles in the current are
   *                  allowed to have.
   * @return Whether the current was sucessfully constructed.
   */
  virtual bool createMode(int icharge, tcPDPtr resonance,
			  FlavourInfo flavour,
			  unsigned int imode,PhaseSpaceModePtr mode,
			  unsigned int iloc,int ires,
			  PhaseSpaceChannel phase, Energy upp );

  /**
   * The particles produced by the current. This just returns the leptons.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual tPDVector particles(int icharge, unsigned int imode, int iq, int ia);
  //@}

  /**
   * Hadronic current. This method is purely virtual and must be implemented in
   * all classes inheriting from this one.
   * @param resonance If specified only include terms with this particle
   * @param flavour Information on the required flavours of the quarks
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVectorE> 
  current(tcPDPtr resonance,
	  FlavourInfo flavour,
	  const int imode, const int ichan,Energy & scale,
	  const tPDVector & outgoing,
	  const vector<Lorentz5Momentum> & momenta,
	  DecayIntegrator::MEOption meopt) const;

  /**
   *   Construct the SpinInfo for the decay products
   */
  void constructSpinInfo(ParticleVector decay) const;

  /**
   * Accept the decay. Checks that this is one of the allowed modes.
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Returns the decay mode number for a given set of particles in the current. 
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

public:

  /**
   * Standard Init function used to initialize the interfaces.
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
   * Private and non-existent assignment operator.
   */
  LeptonNeutrinoCurrent & operator=(const LeptonNeutrinoCurrent &) = delete;

private:
  
  /**
   *  Spinors for the decay products
   */
  mutable vector<Helicity::LorentzSpinor   <SqrtEnergy> > wave_;

  /**
   *  barred spinors for the decay products
   */
  mutable vector<Helicity::LorentzSpinorBar<SqrtEnergy> > wavebar_;
  
};

}


#endif /* HERWIG_LeptonNeutrinoCurrent_H */
