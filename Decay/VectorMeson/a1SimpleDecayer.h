// -*- C++ -*-
//
// a1SimpleDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_a1SimpleDecayer_H
#define HERWIG_a1SimpleDecayer_H
//
// This is the declaration of the a1SimpleDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig/Decay/ResonanceHelpers.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The a1SimpleDecayer class provides a simple model of the decay of the
 * \f$a_1\f$ meson to three pions including \f$\rho\f$ meson intermediate states.
 *
 * @see \ref a1SimpleDecayerInterfaces "The interfaces"
 * defined for a1SimpleDecayer.
 */
class a1SimpleDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  a1SimpleDecayer();

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const tPDVector & outgoing,
	     const vector<Lorentz5Momentum> & momenta,
	     MEOption meopt) const;

  /**
   *   Construct the SpinInfos for the particles produced in the decay
   */
  virtual void constructSpinInfo(const Particle & part,
				 ParticleVector outgoing) const;

  /**
   * Method to return an object to calculate the 3 body partial width.
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;

  /**
   * The matrix element to be integrated for the three-body decays as a function
   * of the invariant masses of pairs of the outgoing particles.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element
   */
  virtual double threeBodyMatrixElement(const int imode , const Energy2 q2,
					const Energy2 s3, const Energy2 s2,
					const Energy2 s1, const Energy  m1,
					const Energy  m2, const Energy m3) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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
  virtual IBPtr clone() const  {return new_ptr(*this);}

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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:
  
  /**
   * The \f$\rho\f$ form factors
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   * @return The form factor
   */
  Complex rhoFormFactor(Energy2 q2,int ires) const {
    Complex output(0.),norm(0.);
    for(unsigned int ix=0;ix<3;++ix) norm += _rhowgts[ix];
    if(ires<0) {
      for(unsigned int ix=0;ix<3;++ix)
	output+=_rhowgts[ix]*Resonance::BreitWignerPWave(q2,_rhomass[ix],_rhowidth[ix],_mpi,_mpi);
    }
    else {
      assert(ires<3);
      output=_rhowgts[ires]*Resonance::BreitWignerPWave(q2,_rhomass[ires],_rhowidth[ires],_mpi,_mpi);
    }
    return output/norm;
  }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  a1SimpleDecayer & operator=(const a1SimpleDecayer &) = delete;

private:

  /**
   * The \f$\rho\f$ masses
   */
  vector<Energy> _rhomass;
  
  /**
   * The \f$\rho\f$ widths
   */
  vector<Energy> _rhowidth;

  /**
   * Weights for the different \f$\rho\f$ resonances 
   */
  vector<double> _rhowgts;

  /**
   *  Use local values of the parameters
   */
  bool _localparameters;

  /**
   *  The overall coupling for the decay
   */
  InvEnergy _coupling;

  /**
   * Maximum weight for the one charged pion channel.
   */
  double _onemax;

  /**
   * Maximum weight for the two charged pion channel.
   */
  double _twomax;

  /**
   * Maximum weight for the three charged pion channel.
   */
  double _threemax;
  
  /**
   * Weights for the channels for the one charged pion channel.
   */
  vector<double> _onewgts;
  
  /**
   * Weights for the channels for the two charged pion channel.
   */
  vector<double> _twowgts;
  
  /**
   * Weights for the channels for the three charged pion channel.
   */
  vector<double> _threewgts;

  /**
   * The pion mass
   */
  Energy _mpi;

  /**
   *  Spin Density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization vectors
   */
  mutable vector<Helicity::LorentzPolarizationVector> _vectors;

};

}

#endif /* HERWIG_a1SimpleDecayer_H */
