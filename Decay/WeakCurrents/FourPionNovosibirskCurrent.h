// -*- C++ -*-
//
// FourPionNovosibirskCurrent.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_FourPionNovosibirskCurrent_H
#define HERWIG_FourPionNovosibirskCurrent_H
//
// This is the declaration of the FourPionNovosibirskCurrent class.
//
#include "WeakDecayCurrent.h"
#include "Herwig/Utilities/Interpolator.h"
#include "Herwig/Utilities/Kinematics.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>FourPionNovosibirskCurrent</code> class implements the decay of the weak 
 * current to 4 pions using the hadronic currents of
 * Comput. Phys. Commun. 146: 139-153, 2002,
 * which is a model based on the \f$e^+e^-\to4\pi\f$ data from Novosibirsk.
 *
 * It should be noted that there were a large number of mistakes in this paper which 
 * were corrected in hep-ph/0312240.
 *
 * @see WeakDecayCurrent
 * @see FourPionDefaultMatrixElement
 * 
 * \author Peter Richardson
 * 
 */
class FourPionNovosibirskCurrent: public WeakDecayCurrent {

  /**
   * The FourPionDefaultMatrixElement class is a friend so it can perform the
   * integration.
   */
  friend class FourPionDefaultMatrixElement;

public:

  /**
   * Default constructor
   */
  FourPionNovosibirskCurrent();

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
   * This version constructs the four pion current.
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
   * The particles produced by the current. This returns the four pions for the
   * current.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual tPDVector particles(int icharge, unsigned int imode, int iq, int ia);
  //@}

  /**
   * Hadronic current. This version calculates the four pion current described above.
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @param meopt Option for the calculation of the matrix element
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVectorE> 
  current(const int imode, const int ichan,Energy & scale, 
	  const ParticleVector & decay, DecayIntegrator::MEOption meopt) const;

  /**
   * Accept the decay. Checks this is one of the four pion modes.
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Return the decay mode number for a given set of particles in the current. 
   * Works out which four pion mode this is.
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

  /**
   * The matrix element to evaluate the \f$a_1\f$ running width.
   * @param iopt The mode
   * @param q2 The mass of the decaying off-shell \f$a_1\f$, \f$q^2\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element squared summed over spins.
   */
  double threeBodyMatrixElement(const int iopt, const Energy2 q2,
				const Energy2 s3, const Energy2 s2, 
				const Energy2 s1, const Energy  m1,
				const Energy  m2, const Energy  m3) const;
  
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
   * Initialize this object to the begining of the run phase.
   */
  virtual void doinitrun();

  /**
   * Check sanity of the object during the setup phase.
   */
  virtual void doupdate();
  //@}

private:

  /**
   * Private and non-existent assignment operator.
   */
  FourPionNovosibirskCurrent & operator=(const FourPionNovosibirskCurrent &) = delete;
      
protected:
  

  /**
   * Initialize the \f$a_1\f$ width.
   * @param iopt Initialization option
   *  (-1 is full initialization and 0 sets up the interpolator for the running width)
   */
  void inita1width(int iopt);
  
  /**
   * Form foactor for the \f$a_1\f$ vertex.
   * @param q2 The scale \f$q^2\f$.
   * @return The \f$a_1\f$ form factor.
   */
  double a1FormFactor(Energy2 q2) const {
    return sqr((1.+_a1massolam2)/(1.+q2*_onedlam2));   
  }

  /**
   * Breit-Wigner for the \f$\sigma\f$ meson
   * @param q2 The scale \f$q^2\f$.
   * @param iopt The pion masses to used (0=\f$\pi^0\f$, 1=\f$\pi^+\f$)
   * @return The Breit-Wigner for the \f$\sigma\f$ meson
   */
  Complex sigmaBreitWigner(Energy2 q2,unsigned int iopt) const;

  /**
   * The \f$a_1\f$ breit wigner.
   * @param q2 The scale \f$q^2\f$.
   * @return The Breit-Wigner for the \f$a_1\f$.
   */
  Complex a1BreitWigner(Energy2 q2) const;

  /**
   * The Breit-Wigner for the \f$\omega\f$.
   * @param q2 The scale \f$q^2\f$.
   * @return The Breit-Wigner for the \f$\omega\f$.
   */
  Complex omegaBreitWigner(Energy2 q2) const;

  /**
   * The Breit-Wigner for the \f$\rho\f$.
   * @param q2 The scale \f$q^2\f$.
   * @return The Breit-Wigner for the \f$\rho\f$.
   */
  Complex rhoBreitWigner(Energy2 q2) const;

  /**
   * Return the \f$a_1\f$ running width.
   * @param q2 The scale \f$q^2\f$.
   * @return The running width.
   */
  Energy a1width(Energy2 q2) const {return (*_a1runinter)(q2);}

  /**
   * The \f$t_1\f$ current used in calculating the current.
   * @param q1 The first momentum.
   * @param q2 The first momentum.
   * @param q3 The first momentum.
   * @param q4 The first momentum.
   * @return The current \f$t_1\f$.
   */
  LorentzVector<complex<Energy5> > 
  t1(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
     Lorentz5Momentum & q3,Lorentz5Momentum & q4) const;

  /**
   * The \f$t_2\f$ current used in calculating the current.
   * @param q1 The first momentum.
   * @param q2 The first momentum.
   * @param q3 The first momentum.
   * @param q4 The first momentum.
   * @param iopt 0 for \f$\sigma\to\pi^+\pi^-\f$ and 1 for \f$\sigma\to\pi^0\pi^0\f$
   * @return The current \f$t_2\f$.
   */
  LorentzVector<complex<Energy5> > 
  t2(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
     Lorentz5Momentum & q3,Lorentz5Momentum & q4,
     unsigned int iopt) const;

  /**
   * The \f$t_3\f$ current used in calculating the current.
   * @param q1 The first momentum.
   * @param q2 The first momentum.
   * @param q3 The first momentum.
   * @param q4 The first momentum.
   * @return The current \f$t_3\f$.
   */
  LorentzVector<complex<Energy5> > 
  t3(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
     Lorentz5Momentum & q3,Lorentz5Momentum & q4) const;

  /**
   * The G functions of hep-ph/0201149
   * @param q2 The scale \f$q^2\f$.
   * @param ichan Which of the four pion channels this is for.
   * @return The G function.
   */
  InvEnergy6 gFunction(Energy2 q2, int ichan) const;

  /**
   * The d parameter in \f$\rho\f$ the propagator.
   */
  Energy2 DParameter() const;

  /**
   * The \f$\frac{dh}{dq^2}\f$ function in the rho propagator evaluated at \f$q^2=m^2\f$.
   */
  double dhdq2Parameter() const;

  /**
   * The h function in the \f$\rho\f$ propagator.
   * @param q The scale.
   * @return The h function.
   */
  Energy2 hFunction(const Energy q) const;

private:
  
  /**
   * Interpolating functions for the G functions of hep-ph/0201149
   */
  //@{
  /**
   * The interpolator for the \f$\omega\f$ current.
   */
  Interpolator<double,Energy>::Ptr _Fomega;

  /**
   * The interpolator for the three charged pion \f$a_1\f$ current. 
   */
  Interpolator<double,Energy>::Ptr _Fthreec;

  /**
   * The interpolator for the one   charged pion \f$a_1\f$ current.
   */
  Interpolator<double,Energy>::Ptr _Fonec;

  /**
   * The interpolator for the \f$\sigma\f$ current.
   */
  Interpolator<double,Energy2>::Ptr _Fsigma;
  //@}

  /**
   * The charged pion mass
   */
  Energy _mpic;

  /**
   * The neutral pion mass
   */
  Energy _mpi0;

  /**
   * The mass of the \f$\rho\f$ for the current.
   */
  Energy _rhomass;

  /**
   * The mass of the \f$a_1\f$ for the current.
   */
  Energy _a1mass;

  /**
   * The mass of the \f$\omega\f$ for the current.
   */
  Energy _omegamass;

  /**
   * The mass of the \f$\sigma\f$ for the current.
   */
  Energy _sigmamass;

  /**
   * The width for the \f$\rho\f$.
   */
  Energy _rhowidth;

  /**
   *  The \f$a_1\f$ width
   */
  Energy _a1width;

  /**
   *  The \f$\omega\f$ width.
   */
  Energy _omegawidth;

  /**
   *  The \f$\sigma\f$ width.
   */
  Energy _sigmawidth;

  /**
   * Mass for the intermediate in the phase-space, this is a technical parameter to
   * improve the phase-space integration efficiency.
   */
  Energy _intmass;

  /**
   * Width for the intermediate in the phase-space, this is a technical parameter to
   * improve the phase-space integration efficiency.
   */
  Energy _intwidth; 

  /**
   * The \f$z\f$ \f$\sigma\f$ coupling.
   */
  Complex _zsigma;

  /**
   * The magnitude of the \f$z\f$ \f$\sigma\f$ coupling.
   */
  double _zmag;

  /**
   * The phase of the \f$z\f$ \f$\sigma\f$ coupling.
   */
  double _zphase;

  /**
   * The mass parameter for the \f$a_1\f$ form-factor.
   */
  Energy2 _lambda2;

  /**
   * The inverse of the mass parameter for the \f$a_1\f$ form-factor.
   */
  InvEnergy2 _onedlam2;

  /**
   *  The physical \f$a_1\f$ mass divided by the mass parameter in the 
   *  \f$a_1\f$ form-factor.
   */
  double _a1massolam2;

  /**
   * The momentum of the  pions in on-shell \f$\sigma\f$ decay which is used
   * in the calculation of the running \f$\sigma\f$ width.
   */
  vector<Energy> _psigma;

  /**
   *  The charged pion mass squared.
   */
  Energy2 _mpic2;

  /**
   * The neutral pion mass squared
   */
  Energy2 _mpi02;

  /**
   *  The h function evaluated at the \f$\rho\f$ mass.
   */
  Energy2 _hm2;

  /**
   * The d parameter for the \f$\rho\f$ width.
   */
  Energy2 _rhoD;

  /**
   * The momentum of the pions produced in on-shell \f$rho\f$ decay.
   */
  Energy _prho;

  /**
   * \f$\frac{dh}{dq^2}\f$ evaluates at \f$q^2=m^2\f$ for the \f$\rho\f$.
   */
  double _dhdq2m2;

  /**
   * Magic number for the omega current.
   */
  InvEnergy _aomega;

  /**
   * Magic number for the three charged pion current.
   */
  InvEnergy _athreec;

  /**
   * Magic number for the one charged pion current
   */
  InvEnergy _aonec;

  /**
   * Magic number for the omega current.
   */
  double _bomega;

  /**
   * Magic number for the three charged pion current.
   */
  double _bthreec;

  /**
   * Magic number for the one charged pion current
   */
  double _bonec;

  /**
   * Magic number for the omega current.
   */
  double _comega;

  /**
   * Magic number for the three charged pion current.
   */
  double _cthreec;

  /**
   * Magic number for the one charged pion current
   */
  double _conec;

  /**
   * magic numbers for the running omega width
   */
  vector<double> _omegaparam;

  /**
   * whether or not to initialize the calculation of the \f$a_1\f$ width
   */
  bool _initializea1;

  /**
   * use local values of the particle masses
   */
  bool _localparameters;

  /**
   * The widths for the interpolation table for the running \f$a_1\f$ width.
   */
  vector<Energy> _a1runwidth;

  /**
   * The \f$q^2\f$ values for the interpolation table for the running \f$a_1\f$ width.
   */
  vector<Energy2> _a1runq2;

  /**
   * The interpolator for the running \f$a_1\f$ width.
   */
  Interpolator<Energy,Energy2>::Ptr _a1runinter;

  /**
   *  The maximum mass of the hadronic system
   */
  Energy _maxmass;

  /**
   *  The maximum mass when the running width was calculated
   */
  Energy _maxcalc;
};

}

#endif /* HERWIG_FourPionNovosibirskCurrent_H */

