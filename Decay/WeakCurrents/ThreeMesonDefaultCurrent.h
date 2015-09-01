// -*- C++ -*-
//
// ThreeMesonDefaultCurrent.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ThreeMesonDefaultCurrent_H
#define HERWIG_ThreeMesonDefaultCurrent_H
//
// This is the declaration of the ThreeMesonDefaultCurrent class.
//
#include "ThreeMesonCurrentBase.h"
#include "Herwig/Utilities/Interpolator.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The ThreeMesonDefaultCurrent class implements the currents from Z.Phys.C58:445 (1992),
 * this paper uses the form from Z.Phys.C48:445 (1990) for the \f$a_1\f$ width and
 * is the default model in TAUOLA.
 *
 *  The following three meson modes are implemented.
 *
 * - \f$    \pi^-  \pi^-    \pi^+ \f$, (imode=0)
 * - \f$    \pi^0  \pi^0    \pi^- \f$, (imode=1)
 * - \f$    K^-   \pi^-    K^+ \f$, (imode=2)
 * - \f$    K^0   \pi^-    \bar{K}^0\f$, (imode=3)
 * - \f$    K^-   \pi^0    K^0 \f$, (imode=4)
 * - \f$    \pi^0  \pi^0    K^- \f$, (imode=5)
 * - \f$    K^-   \pi^-    \pi^+ \f$, (imode=6)
 * - \f$    \pi^-  \bar{K}^0  \pi^0 \f$, (imode=7)
 * - \f$    \pi^-  \pi^0    \eta \f$, (imode=8)
 *
 *  using the currents from TAUOLA
 *
 *
 * @see ThreeMesonCurrentBase,
 * @see WeakDecayCurrent.
 * @see Defaulta1MatrixElement
 * 
 */
class ThreeMesonDefaultCurrent: public ThreeMesonCurrentBase {

  /**
   * The matrix element for the running \f$a_1\f$ width is a friend to 
   * keep some members private.
   */
  friend class Defaulta1MatrixElement;

public:

  /**
   * Default constructor
   */
  ThreeMesonDefaultCurrent();

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
   * This version addes the mesons for the current
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
  //@}

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
  ThreeMesonDefaultCurrent & operator=(const ThreeMesonDefaultCurrent &);

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
   * \f$a_1\f$ Breit-Wigner
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @return The Breit-Wigner
   */
  Complex a1BreitWigner(Energy2 q2) const  {
    Complex ii(0.,1.);
    Energy2 m2(_a1mass*_a1mass);
    Energy  q(sqrt(q2));
    return m2/(m2-q2-ii*q*a1Width(q2));
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
   * The \f$a_1\f$ running width
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @return The \f$a_1\f$ running width.
   */
  Energy a1Width(Energy2 q2) const {
    Energy output;
    if(!_a1opt) output = _a1mass*_a1width*g(q2)/g(sqr(_a1mass))/sqrt(q2);
    else        output = (*_a1runinter)(q2);
    return output;
  }
  
  /**
   *  The \f$g(Q^2)\f$ function of Kuhn and Santamaria
   */
  double g(Energy2 q2) const {
    double output;
    if(q2 < 9.*sqr(_mpi)) {
      output=0.;
    }
    else if(q2 < sqr(_rhoF123masses[0]+_mpi)) {
      double diff = (q2-9.*sqr(_mpi))/GeV2;
      
      output = 4.1*sqr(diff)*diff*(1.-3.3*diff+5.8*sqr(diff));
    }
    else {
      double ratio = q2/GeV2;
      output = ratio*(1.623+10.38/ratio-9.32/sqr(ratio)+0.65/(ratio*sqr(ratio)));
    }
    return output;
  }

  /**
   * Initialize the \f$a_1\f$ running width
   * @param iopt Initialization option (-1 full calculation, 0 set up the interpolation)
   */
  void inita1Width(int iopt);

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

  /**
   * Initialize the running \f$a_1\f$ width.
   */
  bool _initializea1;
  
  /**
   * The mass of the \f$a_1\f$ resonances.
   */
  Energy _a1mass;

  /**
   * The width of the \f$a_1\f$ resonances.
   */
  Energy _a1width;

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
   * use local values of the \f$\rho\f$ masses and widths
   */
  bool _rhoparameters;

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
   * use local values of the \f$K^*\f$ resonances masses and widths
   */
  bool _kstarparameters;

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
  
  /**
   * Use local values of the \f$a_1\f$ parameters
   */
  bool _a1parameters;
  
  /**
   * Use local values of the \f$K_1\f$ parameters
   */
  bool _k1parameters;

  /**
   * Option for the \f$a_1\f$ width
   */
  bool _a1opt;

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

#endif /* THEPEG_ThreeMesonDefaultCurrent_H */
