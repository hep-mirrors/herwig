// -*- C++ -*-
//
// OneKaonTwoPionDefaultCurrent.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_OneKaonTwoPionDefaultCurrent_H
#define HERWIG_OneKaonTwoPionDefaultCurrent_H
//
// This is the declaration of the OneKaonTwoPionDefaultCurrent class.
//
#include "WeakCurrent.h"
#include "Herwig/Utilities/Interpolator.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/ResonanceHelpers.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The OneKaonTwoPionDefaultCurrent class implements the currents from Z.Phys.C58:445 (1992),
 * this paper uses the form from Z.Phys.C48:445 (1990) for the \f$a_1\f$ width and
 * is the default model in TAUOLA.
 *
 *  The following three meson modes are implemented.
 *
 * - \f$    \pi^0  \pi^0    K^- \f$, (imode=5)
 * - \f$    K^-   \pi^-    \pi^+ \f$, (imode=6)
 * - \f$    \pi^-  \bar{K}^0  \pi^0 \f$, (imode=7)
 *
 *  using the currents from TAUOLA
 *
 *
 * @see WeakCurrent
 * @see Defaulta1MatrixElement
 * 
 */
class OneKaonTwoPionDefaultCurrent: public WeakCurrent {

public:

  /**
   * Default constructor
   */
  OneKaonTwoPionDefaultCurrent();

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
	  IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
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
			  IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
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

  /**
   * Helper class for form factors
   */
  struct FormFactors {

    /**
     * @param F1 The \f$F_1\f$ form factor
     */
    complex<InvEnergy>  F1;
    
    /**
     * @param F2 The \f$F_2\f$ form factor
     */
    complex<InvEnergy>  F2;
    
    /**
     * @param F3 The \f$F_3\f$ form factor
     */
    complex<InvEnergy>  F3; 
    
    /**
     * @param F4 The \f$F_4\f$ form factor
     */
    complex<InvEnergy>  F4;
    
    /**
     * @param F5 The \f$F_5\f$ form factor
     */
    complex<InvEnergy3> F5;

    /**
     *  Constructor
     * @param f1 The \f$F_1\f$ form factor
     * @param f2 The \f$F_2\f$ form factor
     * @param f3 The \f$F_3\f$ form factor
     * @param f4 The \f$F_4\f$ form factor
     * @param f5 The \f$F_5\f$ form factor
     */    
    FormFactors(complex<InvEnergy>  f1 = InvEnergy(), 
		complex<InvEnergy>  f2 = InvEnergy(),
		complex<InvEnergy>  f3 = InvEnergy(),
		complex<InvEnergy>  f4 = InvEnergy(),
		complex<InvEnergy3> f5 = InvEnergy3())
      : F1(f1), F2(f2), F3(f3), F4(f4), F5(f5) {}
  };
  
protected:

  /**
   * Can the current handle a particular set of mesons. 
   * As this current includes all the allowed modes this is always true.
   */
  virtual bool acceptMode(int) const;

  /**
   * Calculate the form factor for the current. Implements the form factors
   * described above.
   * @param ichan The phase space channel
   * @param imode The mode
   * @param q2 The scale \f$q^2\f$ for the current.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   */
  virtual FormFactors calculateFormFactors(const int ichan, const int imode,
					   Energy2 q2,
					   Energy2 s1, Energy2 s2, Energy2 s3) const;

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
  OneKaonTwoPionDefaultCurrent & operator=(const OneKaonTwoPionDefaultCurrent &);

private:
  
  /**
   * The \f$\rho\f$ Breit-Wigner for the \f$F_{1,2,3}\f$ form factors.
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   * @return The Breit-Wigner 
   */
  Complex BrhoF123(Energy2 q2,int ires) const {
    Complex output(0.),norm(0.);
    for(unsigned int ix=0,N=min(3,int(_rhoF123wgts.size()));ix<N;++ix) {
      norm+=_rhoF123wgts[ix];
    }
    if(ires<0) {
      for(unsigned int ix=0,N=min(3,int(_rhoF123wgts.size()));ix<N;++ix) {
	output+=_rhoF123wgts[ix]*rhoKBreitWigner(q2,0,ix);
      }
    }
    else {
      unsigned int temp(ires);
      if(temp<_rhoF123wgts.size()&&temp<3)
	output=_rhoF123wgts[temp]*rhoKBreitWigner(q2,0,temp);
      else
	output=0.;
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
    Complex output(0.),norm(0.);
    for(unsigned int ix=0,N=min(3,int(_rhoF5wgts.size()));ix<N;++ix) {
      norm+=_rhoF5wgts[ix];
    }
    if(ires<0) {
      for(unsigned int ix=0,N=min(3,int(_rhoF5wgts.size()));ix<N;++ix) {
	output+=_rhoF5wgts[ix]*rhoKBreitWigner(q2,1,ix);
      }
    }
    else {
      unsigned int temp(ires);
      if(temp<_rhoF5wgts.size()&&temp<3) {
	output=_rhoF5wgts[temp]*rhoKBreitWigner(q2,1,temp);
      }
    }
    return output/norm;
  }

  /**
   * The \f$K^*\f$ Breit-Wigner for the \f$F_{1,2,3}\f$ form factors.
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   * @return The Breit-Wigner 
   */
  Complex BKstarF123(Energy2 q2,int ires) const {
    Complex output(0.),norm(0.);
    for(unsigned int ix=0,N=min(3,int(_kstarF123wgts.size()));ix<N;++ix) {
      norm+=_kstarF123wgts[ix];
    }
    if(ires<0) {
      for(unsigned int ix=0,N=min(3,int(_kstarF123wgts.size()));ix<N;++ix) {
	output+=_kstarF123wgts[ix]*rhoKBreitWigner(q2,2,ix);
      }
    }
    else {
      unsigned int temp(ires);
      if(temp<_kstarF123wgts.size()&&temp<3) {
	output=_kstarF123wgts[temp]*rhoKBreitWigner(q2,2,temp);
      }
    }
    return output/norm;
  }

  /**
   * The \f$K^*\f$ Breit-Wigner for the \f$F_5\f$ form factors.
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   * @return The Breit-Wigner 
   */
  Complex BKstarF5(Energy2 q2,int ires) const {
    Complex output(0.),norm(0.);
    for(unsigned int ix=0,N=min(3,int(_kstarF5wgts.size()));ix<N;++ix) {
      norm+=_kstarF5wgts[ix];
    }
    if(ires<0) {
      for(unsigned int ix=0,N=min(3,int(_kstarF5wgts.size()));ix<N;++ix) {
	output+=_kstarF5wgts[ix]*rhoKBreitWigner(q2,3,ix);
      }
    }
    else {
      unsigned int temp(ires);
      if(temp<_kstarF5wgts.size()&&temp<3) {
	output=_kstarF5wgts[ires]*rhoKBreitWigner(q2,3,temp);
      }
    }
    return output/norm;
  }
  
  /**
   * Mixed Breit Wigner for the \f$F_5\f$ form factor
   * @param si The scale \f$s_1\f$.
   * @param sj The scale \f$s_2\f$.
   * @param ires Which resonances to use
   * @return The mixed Breit-Wigner
   */
  Complex FKrho(Energy2 si,Energy2 sj,int ires) const {
    Complex output;
    if(ires<0){output = _rhoKstarwgt*BKstarF123(si,-1)+BrhoF123(sj,-1);}
    else if(ires%2==0){output= _rhoKstarwgt*BKstarF123(si,ires/2);}
    else if(ires%2==1){output=BrhoF123(sj,ires/2);}
    output /=(1.+_rhoKstarwgt);
    return output;
  }
  
  /**
   * The \f$K_1\f$ Breit-Wigner
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @return The Breit-Wigner
   */
  Complex K1BreitWigner(Energy2 q2) const {
    Energy2 m2 = sqr(_k1mass);
    Complex ii(0.,1.);
    complex<Energy2> fact(m2 - ii*_k1mass*_k1width);
    return fact/(fact-q2);
  }

  /**
   * Breit-Wigners for the \f$\rho\f$ and \f$K^*\f$.
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner.
   * @param itype The type of Breit-Wigner, \e i.e. which masses and widths to use.x
   * @param ires Which multiplet to use.
   */
  Complex rhoKBreitWigner(Energy2 q2,unsigned int itype,unsigned int ires) const;

private:
  
  /**
   * Parameters for the \f$\rho\f$ Breit-Wigner in the
   * \f$F_{1,2,3}\f$ form factors.
   */
  vector<double> _rhoF123wgts;

  /**
   * Parameters for the \f$K^*\f$ Breit-Wigner in the
   * \f$F_{1,2,3}\f$ form factors.
   */
  vector<double> _kstarF123wgts;
  
  /**
   * Parameters for the \f$\rho\f$ Breit-Wigner in the
   * \f$F_5\f$ form factors.
   */
  vector<double> _rhoF5wgts;

  /**
   * Parameters for the \f$K^*\f$ Breit-Wigner in the
   * \f$F_5\f$ form factors.
   */
  vector<double> _kstarF5wgts;
  
  /**
   * The relative weight of the \f$\rho\f$ and \f$K^*\f$ where needed.
   */
  double _rhoKstarwgt;

  /**
   * The mass of the \f$aK1\f$ resonances.
   */
  Energy _k1mass;

  /**
   * The width of the \f$K_1\f$ resonances.
   */
  Energy _k1width;

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

  /**
   * The \f$K^*\f$ masses for the \f$F_{1,2,3}\f$ form factors.
   */
  vector<Energy> _kstarF123masses;

  /**
   * The \f$K^*\f$ masses for the \f$F_5\f$ form factors.
   */
  vector<Energy> _kstarF5masses;

  /**
   * The \f$K^*\f$ widths for the \f$F_{1,2,3}\f$ form factors.
   */
  vector<Energy> _kstarF123widths;

  /**
   * The \f$K^*\f$ widths for the \f$F_5\f$ form factors.
   */
  vector<Energy> _kstarF5widths;
};

}

#endif /* HERWIG_OneKaonTwoPionDefaultCurrent_H */
