// -*- C++ -*-
//
// TwoKaonOnePionCurrent.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TwoKaonOnePionCurrent_H
#define HERWIG_TwoKaonOnePionCurrent_H
//
// This is the declaration of the TwoKaonOnePionCurrent class.
//

#include "WeakCurrent.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The TwoKaonOnePionCurrent class implements the model of M. Finkemeier 
 * and E.~Mirkes, Z. Phys. C 69 (1996) 243 [arXiv:hep-ph/9503474],
 * for the weak current for three mesons where at least one of the mesons is
 * a kaon.
 *
 * \ingroup Decay
 *
 *  This is the base class for the three meson decays of the weak current.
 *  It is designed so that the currents for the following modes can be implemented
 *  in classes inheriting from this
 * - \f$    K^-   \pi^-    K^+ \f$, (imode=0)
 * - \f$    K^0   \pi^-    \bar{K}^0\f$, (imode=1)
 * - \f$    K^-   \pi^0    K^0 \f$, (imode=2)
 * - \f$    \pi^0  \pi^0    K^- \f$, (imode=3)
 * - \f$    K^-   \pi^-    \pi^+ \f$, (imode=4)
 * - \f$    \pi^-  \bar{K}^0  \pi^0 \f$, (imode=5)
 * - \f$    \pi^-  \pi^0    \eta \f$, (imode=6)
 *
 * obviously there are other modes with three pseudoscalar mesons for the decay
 * of the weak current but this model original came from \f$\tau\f$ decay where
 * these are the only modes. However one case which is important is the inclusion
 * of the mixing in the neutral kaon sector for which we include the additional
 * currents
 * - \f$    K^0_S \pi^- K^0_S\f$, (imode=9)
 * - \f$    K^0_L \pi^- K^0_L\f$, (imode=10)
 * - \f$    K^0_S \pi^- K^0_L\f$, (imode=11)
 *
 *  In this case the current is given by
 *  \f[ J^\mu = \left(g^{\mu\nu}-\frac{q^\mu q^\nu}{q^2}\right)
 *   \left[F_1(p_2-p_3)^\mu +F_2(p_3-p_1)^\mu+F_3(p_1-p_2)^\mu\right]
 *  +q^\mu F_4
 *  +F_5\epsilon^{\mu\alpha\beta\gamma}p_1^\alpha p_2^\beta p_3^\gamma
 *  \f]
 * where
 * - \f$p_{1,2,3}\f$ are the momenta of the mesons in the order given above.
 * - \f$F_1,F_2,F_3,F_4,F_5\f$ are the form factors which must be 
 *  calculated in the calculateFormFactors member which should be implemented
 * in classes inheriting from this.
 *
 * @see WeakCurrent.
 *  
 * \author Peter Richardson
 * @see \ref TwoKaonOnePionCurrentInterfaces "The interfaces"
 * defined for TwoKaonOnePionCurrent.
 */
class TwoKaonOnePionCurrent: public WeakCurrent {

public:

  /**
   * The default constructor.
   */
  TwoKaonOnePionCurrent();

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
  
  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;
  
  /**
   * the matrix element for the \f$a_1\f$ decay to calculate the running width
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The mass of the decaying off-shell \f$a_1\f$, \f$q^2\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element squared summed over spins.
   */
  double threeBodyMatrixElement(const int imode,  const Energy2 q2,
				const Energy2 s3, const Energy2 s2, 
				const Energy2 s1, const Energy  m1, 
				const Energy  m2, const Energy  m3) const;
  
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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TwoKaonOnePionCurrent & operator=(const TwoKaonOnePionCurrent &) = delete;

private:

  /**
   *  The \f$\rho\f$ lineshape for the axial-vector terms
   * @param q2 The scale \f$q^2\f$ for the lineshape
   * @param ires Which \f$\rho\f$ multiplet
   */
  Complex Trho1(Energy2 q2,int ires) const {
    if(ires>=int(_rho1wgts.size())) return 0.;
    double norm = std::accumulate(_rho1wgts.begin(),_rho1wgts.end(),0.);
    unsigned int imin=0,imax=_rho1wgts.size();
    if(ires>0) {
      imin=ires;
      imax=imin+1;
    }
    Complex output(0.);
    for(unsigned int ix=imin;ix<imax;++ix)
      output+=_rho1wgts[ix]*
	Resonance::BreitWignerPWave(q2,_rho1mass[ix],_rho1width[ix],_mpi,_mpi);
    return output/norm;
  }

  /**
   *  The \f$\rho\f$ lineshape for the vector terms
   * @param q2 The scale \f$q^2\f$ for the lineshape
   * @param ires Which \f$\rho\f$ multiplet
   */
  Complex Trho2(Energy2 q2,int ires) const {
    if(ires>=int(_rho2wgts.size())) return 0.;
    double norm = std::accumulate(_rho2wgts.begin(),_rho2wgts.end(),0.);
    unsigned int imin=0,imax=_rho2wgts.size();
    if(ires>0) {
      imin=ires;
      imax=imin+1;
    }
    Complex output(0.);
    for(unsigned int ix=imin;ix<imax;++ix)
      output+=_rho2wgts[ix]*
	Resonance::BreitWignerPWave(q2,_rho2mass[ix],_rho2width[ix],_mpi,_mpi);
    return output/norm;
  }

  /**
   *  The \f$K^*\f$ lineshape for the axial-vector terms
   * @param q2 The scale \f$q^2\f$ for the lineshape
   * @param ires Which \f$K^*\f$ multiplet
   */
  Complex TKstar1(Energy2 q2,int ires) const {
    if(ires>=int(_kstar1wgts.size())) return 0.;
    double norm = std::accumulate(_kstar1wgts.begin(),_kstar1wgts.end(),0.);
    unsigned int imin=0,imax=_kstar1wgts.size();
    if(ires>0) {
      imin=ires;
      imax=imin+1;
    }
    Complex output(0.);
    for(unsigned int ix=imin;ix<imax;++ix)
      output+=_kstar1wgts[ix]*
	Resonance::BreitWignerPWave(q2,_kstar1mass[ix],_kstar1width[ix],_mK,_mpi);
    return output/norm;
  }

  /**
   * \f$a_1\f$ Breit-Wigner
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @return The Breit-Wigner
   */
  Complex a1BreitWigner(Energy2 q2) const {
    Complex ii(0.,1.);
    Energy2 m2(_a1mass*_a1mass);
    Energy  q(sqrt(q2));
    Energy gamma = !_a1opt ?
      _a1mass*_a1width*Resonance::ga1(q2)/Resonance::ga1(_a1mass*_a1mass)/sqrt(q2) : (*_a1runinter)(q2);
    return m2/(m2-q2-ii*q*gamma);
  }

  /**
   * Initialize the \f$a_1\f$ running width
   * @param iopt Initialization option (-1 full calculation, 0 set up the interpolation)
   */
  void inita1Width(int iopt);

  /**
   *  The \f$T_\omega\f$ function
   * @param q2 The scale
   * @param ires the resonance
   */
  Complex Tomega(Energy2 q2, int ires) const;

  /**
   *  The \f$\omega\f$ and \f$\phi\f$ Breit-Wigner
   * @param q2 The scale
   * @param ires the resonance
   */
  Complex OmegaPhiBreitWigner(Energy2 q2, unsigned int ires) const {
    Energy2 m2,mg;
    if(ires==0) {
      m2=sqr(_omegamass);
      mg=_omegamass*_omegawidth;
    }
    else {
      m2=sqr(_phimass);
      mg=_phimass*_phiwidth;
    }
    return (-m2+Complex(0.,1.)*mg)/(q2-m2+Complex(0.,1.)*mg);
  }

  /**
   * The \f$\omega-\phi\f$ \f$K^*\f$ form-factor for the \f$F_5\f$ form-factor
   * @param s1 The scale \f$s_1\f$.
   * @param s2 The scale \f$s_2\f$.
   * @param ires Which resonances to use
   * @return The mixed Breit-Wigner
   */
  Complex TOmegaKStar(Energy2 s1,Energy2 s2,int ires) const;

private:

  /**
   *  Parameters for the \f$\rho\f$ in the axial-vector terms
   */
  //@{
  /**
   *  Weight for the different resonances
   */
  vector<double> _rho1wgts;

  /**
   *  Masses
   */
  vector<Energy> _rho1mass;

  /**
   *  Widths
   */
  vector<Energy> _rho1width;
  //@}

  /**
   *  Parameters for the \f$\rho\f$ in the vector terms
   */
  //@{
  /**
   *  Weight for the different resonances
   */
  vector<double> _rho2wgts;

  /**
   *  Masses
   */
  vector<Energy> _rho2mass;

  /**
   *  Widths
   */
  vector<Energy> _rho2width;
  //@}

  /**
   *  Parameters for the \f$K^*\f$ in the axial-vector terms
   */
  //@{
  /**
   *  Weight for the different resonances
   */
  vector<double> _kstar1wgts;

  /**
   *  Masses
   */
  vector<Energy> _kstar1mass;

  /**
   *  Widths
   */
  vector<Energy> _kstar1width;
  //@}

  /**
   *  Parameters for the three meson resonances
   */
  //@{
  /**
   * The mass of the \f$a_1\f$ resonances.
   */
  Energy _a1mass;
  
  /**
   * The width of the \f$a_1\f$ resonances.
   */
  Energy _a1width;

  /**
   * The \f$a_1\f$ width for the running \f$a_1\f$ width calculation.
   */
  vector<Energy>  _a1runwidth;

  /**
   * The \f$q^2\f$ for the running \f$a_1\f$  width calculation.
   */
  vector<Energy2> _a1runq2;

  /**
   * The interpolator for the running \f$a_1\f$ width calculation.
   */
  Interpolator<Energy,Energy2>::Ptr _a1runinter;
  //@}

  /**
   *  Parameters for the \f$T_\omega\f$ function
   */
  //@{
  /**
   *  Mixing parameter
   */
  double _epsomega;

  /**
   *  Mass of the \f$\omega\f$
   */
  Energy _omegamass;

  /**
   *  Width of the \f$\omega\f$
   */
  Energy _omegawidth;

  /**
   *  Mass of the \f$\phi\f$
   */
  Energy _phimass;

  /**
   *  Width of the \f$\phi\f$
   */
  Energy _phiwidth;
  //@}

  /**
   * The relative weight of the \f$\omega-\phi\f$ and \f$K^*\f$ where needed.
   */
  double _omegaKstarwgt;

  /**
   * The pion decay constant, \f$f_\pi\f$.
   */
  Energy _fpi;

  /**
   * The pion mass
   */
  Energy _mpi;

  /**
   * The kaon mass
   */
  Energy _mK;

  /**
   *  Initialization switches
   */
  //@{
  /**
   * Initialize the running \f$a_1\f$ width.
   */
  bool _initializea1;

  /**
   * Option for the \f$a_1\f$ width
   */
  bool _a1opt;
  //@}

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

#endif /* HERWIG_TwoKaonOnePionCurrent_H */
