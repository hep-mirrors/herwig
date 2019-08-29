// -*- C++ -*-
//
// EtaPiPiDefaultCurrent.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_EtaPiPiDefaultCurrent_H
#define HERWIG_EtaPiPiDefaultCurrent_H
//
// This is the declaration of the EtaPiPiDefaultCurrent class.
//
#include "WeakCurrent.h"
#include "Herwig/Utilities/Interpolator.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/ResonanceHelpers.h"
#include <numeric>

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The EtaPiPiDefaultCurrent class implements the current from Z.Phys.C58:445 (1992),
 * for \f$    \pi^-  \pi^0    \eta \f$.
 *
 *
 * @see WeakCurrent
 * 
 */
class EtaPiPiDefaultCurrent: public WeakCurrent {

public:

  /**
   * Default constructor
   */
  EtaPiPiDefaultCurrent();

  /**
   * Hadronic current. This method is purely virtual and must be implemented in
   * all classes inheriting from this one.
   * @param resonance If specified only include terms with this particle
   * @param Itotal    If specified the total isospin of the current
   * @param I3        If specified the thrid component of isospin
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
	  IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3, Strangeness::Strange S,
	  const int imode, const int ichan,Energy & scale,
	  const tPDVector & outgoing,
	  const vector<Lorentz5Momentum> & momenta,
	  DecayIntegrator::MEOption meopt) const;

  /**
   * Accept the decay. Checks the mesons against the list.
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Return the decay mode number for a given set of particles in the current. 
   * Checks the mesons against the list.
   * @param id The id's of the particles in the current.
   * @return The number of the mode
   */
  virtual unsigned int decayMode(vector<int> id);

  /**
   * The particles produced by the current. This returns the mesons for the mode.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual tPDVector particles(int icharge, unsigned int imode, int iq, int ia);

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
   * Complete the construction of the decay mode for integration.classes inheriting
   * from this one.
   * This method is purely virtual and must be implemented in the classes inheriting
   * from WeakCurrent.
   * @param icharge   The total charge of the outgoing particles in the current.
   * @param resonance If specified only include terms with this particle
   * @param Itotal    If specified the total isospin of the current
   * @param I3        If specified the thrid component of isospin
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
			  IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3, Strangeness::Strange S,
			  unsigned int imode,PhaseSpaceModePtr mode,
			  unsigned int iloc,int ires,
			  PhaseSpaceChannel phase, Energy upp );
  //@}

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

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

private:

  /**
   * Private and non-existent assignment operator.
   */
  EtaPiPiDefaultCurrent & operator=(const EtaPiPiDefaultCurrent &) = delete;

private:
  
  /**
   * The \f$\rho\f$ Breit-Wigner for the \f$F_{1,2,3}\f$ form factors.
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   * @return The Breit-Wigner 
   */
  Complex BrhoF123(Energy2 q2,int ires) const {
    Complex output(0.);
    Complex norm = std::accumulate(_rhoF123wgts.begin(),_rhoF123wgts.end(),Complex(0.));
    if(ires<0) {
      for(unsigned int ix=0;ix<_rhoF123wgts.size();++ix) {
	output+=_rhoF123wgts[ix]*
	  Resonance::BreitWignerPWave(q2,_rhoF123masses[ix],
				      _rhoF123widths[ix],_mpi,_mpi);
      }
    }
    else {
      assert(ires<=int(_rhoF123wgts.size()));
      output=_rhoF123wgts[ires]*
	Resonance::BreitWignerPWave(q2,_rhoF123masses[ires],
				    _rhoF123widths[ires],_mpi,_mpi);
    }
    return output/norm;
  }

  /**
   * The \f$\rho\f$ Breit-Wigner for the \f$F_5\f$ form factors.
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   * @return The Breit-Wigner 
   */
  Complex BrhoF5(Energy2 q2,int ires) const {
    Complex output(0.);
    Complex norm = std::accumulate(_rhoF5wgts.begin(),_rhoF5wgts.end(),Complex(0.0));
    if(ires<0) {
      for(unsigned int ix=0;ix<_rhoF5wgts.size();++ix) {
	output+=_rhoF5wgts[ix]*
	  Resonance::BreitWignerPWave(q2,_rhoF5masses[ix],
				      _rhoF5widths[ix],_mpi,_mpi);
      }
    }
    else {
      assert(ires<=int(_rhoF123wgts.size()));
      output=_rhoF5wgts[ires]*
	Resonance::BreitWignerPWave(q2,_rhoF5masses[ires],
				    _rhoF5widths[ires],_mpi,_mpi);
    }
    return output/norm;
  }

private:
  
  /**
   * Parameters for the \f$\rho\f$ Breit-Wigner in the
   * \f$F_{1,2,3}\f$ form factors.
   */
  vector<double> _rhoF123wgts;
  
  /**
   * Parameters for the \f$\rho\f$ Breit-Wigner in the
   * \f$F_5\f$ form factors.
   */
  vector<double> _rhoF5wgts;

  /**
   * The pion decay constant, \f$f_\pi\f$.
   */
  Energy _fpi;

  /**
   * The pion mass
   */
  Energy _mpi;

  /**
   * The \f$\rho\f$ masses for the \f$F_{1,2,3}\f$ form factors.
   */
  vector<Energy> _rhoF123masses;

  /**
   * The \f$\rho\f$ masses for the \f$F_5\f$ form factors.
   */
  vector<Energy> _rhoF5masses;

  /**
   * The \f$\rho\f$ widths for the \f$F_{1,2,3}\f$ form factors.
   */
  vector<Energy> _rhoF123widths;

  /**
   * The \f$\rho\f$ widths for the \f$F_5\f$ form factors.
   */
  vector<Energy> _rhoF5widths;
};

}

#endif /* HERWIG_EtaPiPiDefaultCurrent_H */
