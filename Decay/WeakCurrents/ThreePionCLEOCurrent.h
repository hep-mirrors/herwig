// -*- C++ -*-
#ifndef THEPEG_ThreePionCLEOCurrent_H
#define THEPEG_ThreePionCLEOCurrent_H
//
// This is the declaration of the ThreePionCLEOCurrent class.
//
#include "ThreeMesonCurrentBase.h"
#include "ThreePionCLEOCurrent.fh"
#include "Herwig++/Utilities/Interpolator.h"
#include "Herwig++/Decay/ThreeBodyIntegrator.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>ThreePionCLEOCurrent</code> class implements
 *  the decay of the weak current to three pions
 *  using the currents from CLEO Phys. Rev. D 61,012002. This is
 *  a model including two \f$\rho\f$ mesons in both \f$s\f$ and \f$p\f$ wave,
 *  a \f$\sigma\f$, the \f$f_2\f$ and \f$f_0(1370)\f$.
 *
 * The form factors for the \f$a_1^+ \to \pi^0 \pi^0 \pi^+\f$ mode are
 *
 * \f[F_1=\sum_k\left\{g^P_{\rho_k}B_{\rho_k}^P(s_1)
 *          -\frac{g^D_{\rho_k}}3B_{\rho_k}^P(s_2)
 *                        \left((s_3-m_{\pi^+}^2)-(s_1-m_{\pi^0}^2)\right)\right\}
 *     +\frac23\left(g_\sigma B^S_\sigma(s_3)+g_{f_0}B^S_{f_0}(s_3)\right)     
 * +\frac{g_{f_2}}{18s_3}(q^2-m_{\pi^+}^2+s_3)(4m_{\pi^0}^2-s_3)B^D_{f_2}(s_3)
 *\f]
 *
 * \f[F_2=\sum_k\left\{-\frac13g^P_{\rho_k}B_{\rho_k}^P(s_2)
 *         -g^D_{\rho_k}B_{\rho_k}^P(s_1)
 *                       \left((s_3-m_{\pi^+}^2)-(s_2-m_{\pi^0}^2)\right)\right\}
 *     +\frac23\left(g_\sigma B^S_\sigma(s_3)+g_{f_0}B^S_{f_0}(s_3)\right)
 * +\frac1{18s_3}g_{f_2}(q^2-m_{\pi^+}^2+s_3)(4m_{\pi^0}^2-s_3)B^D_{f_2}(s_3)
 *\f]
 *
 * \f[F_3=\sum_k g^D_{\rho_k}\left\{ 
 *     -\frac13B_{\rho_k}^P(s_1)\left((s_3-m_{\pi^+}^2)-(s_2-m_{\pi^0}^2)\right)
 *     +\frac13B_{\rho_k}^P(s_2)\left((s_3-m_{\pi^+}^2)-(s_1-m_{\pi^0}^2)\right)\right\}
 *  -\frac{g_{f_2}}2(s_1-s_2)B^D_{f_2}(s_3)\f]
 *  The form factors for  \f$a_1^+\to \pi^+ \pi^+ \pi^-\f$ mode
 *
 * \f[F_1=\sum_k\left\{-g^P_{\rho_k}B_{\rho_k}^P(s_1)
 *                       -\frac{g^D_{\rho_k}}3B_{\rho_k}^P(s_2)(s_1-s_3)\right\}
 *    -\frac23\left(g_\sigma B^S_\sigma(s_2)+g_{f_0} B^S_{f_0}(s_2)\right) 
 *    +g_{f_2}\left(\frac12(s_3-s_2)B^D_{f_2}(s_1)
 *    -\frac1{18s_2}(4m_{\pi^+}^2-s_2)(q^2+s_2-m_{\pi^+}^2)B^D_{f_2}(s_2)\right)\f]
 *
 * \f[F_2=\sum_k\left\{-g^P_{\rho_k}B_{\rho_k}^P(s_2)
 *                       -\frac{g^D_{\rho_k}}3B_{\rho_k}^P(s_1)(s_2-s_3)\right\}
 *    -\frac23\left(g_\sigma B^S_\sigma(s_1)+g_{f_0} B^S_{f_0}(s_1)\right)
 *    +g_{f_2}\left(\frac12(s_3-s_1)B^D_{f_2}(s_2)
 *    -\frac1{18s_1}(4m_{\pi^+}^2-s_1)(q^2+s_1-m_{\pi^+}^2)B^D_{f_2}(s_1)\right)\f]
 *
 * \f[F_3=\sum_k
 *     -g^D_{\rho_k}\left( \frac13(s_2-s_3)B_{\rho_k}^P(s_1)
 *                        -\frac13(s_1-s_3)B_{\rho_k}^P(s_2)\right)
 *   -\frac23\left(g_\sigma B^S_\sigma(s_1)+g_{f_0}B^S_{f_0}(s_1)\right)
 *   +\frac23\left(g_\sigma B^S_\sigma(s_2)+g_{f_0}B^S_{f_0}(s_2)\right)\f]
 *\f[
 *   +g_{f_2}\left(-\frac1{18s_1}(4m_{\pi^+}^2-s_1)(q^2+s_1-m_{\pi^+}^2)B^D_{f_2}(s_1)
 *            +\frac1{18s_2}(4m_{\pi^+}^2-s_2)(q^2+s_2-m_{\pi^+}^2)B^D_{f_2}(s_2)\right)\f]
 *
 * where
 *
 * - \f$g_{f_2}\f$ is the coupling of the \f$f_2\f$ to the \f$a_1\f$
 * - \f$g_{f_0}\f$ is the coupling of the \f$f_0(1370)\f$ to the \f$a_1\f$
 * - \f$g_{\sigma}\f$ is the coupling of the \f$\sigma\f$ to the \f$a_1\f$
 * - \f$g^P_{\rho_k}\f$ is the \f$p\f$-wave coupling of the \f$\rho_k\f$ multiplet
 *     to the \f$a_1\f$.
 * - \f$g^D_{\rho_k}\f$ is the \f$d\f$-wave coupling of the \f$\rho_k\f$ multiplet
 *     to the \f$a_1\f$.
 * - \f$s_3=m^2_{12}\f$ is the invariant mass squared of particles 1 and 2.
 * - \f$s_2=m^2_{13}\f$ is the invariant mass squared of particles 1 and 3.
 * - \f$s_1=m^2_{23}\f$ is the invariant mass squared of particles 2 and 3.
 *
 * The Breit-Wigner factors are given by
    \f$B^L_Y(s_i) = \frac{m^2_Y}{m^2_Y-s_i+im_Y\Gamma^{Y,L}(s_i)}\f$
 * where
 * \f$\Gamma^{Y,L}(s_i) = \Gamma^Y\left(\frac{p(s_i)}{p(M_Y}\right)^{2L+1}\frac{m_Y}{\sqrt{s_i}}\f$
 * \f$m_Y\f$ and \f$\Gamma^Y\f$ are the mass and width of the particle \f$Y\f$ 
 * respectively. \f$p(s_i)\f$ is the momentum of the outgoing pion in the 
 * rest frame of the resonance \f$Y\f$.
 *
 * @see ThreeMesonCurrentBase
 * @see a1ThreePionCLEODecayer
 * @see ThreePionCLEOa1MatrixElement
 * 
 */
class ThreePionCLEOCurrent: public ThreeMesonCurrentBase {

  /**
   * The matrix element for the running \f$a_1\f$ width is a friend to 
   * keep some members private.
   */
  friend class ThreePionCLEOa1MatrixElement;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline ThreePionCLEOCurrent();

  /**
   * Copy constructor
   */
  inline ThreePionCLEOCurrent(const ThreePionCLEOCurrent &);
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

  /**
   * The particles produced by the current. This returns the mesons for the mode.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual PDVector particles(int icharge, unsigned int imode, int iq, int ia);
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
   * the matrix element for the a1 decay to calculate the running width
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
  inline double a1MatrixElement(int iopt,Energy2 q2, Energy2 s3,Energy2 s2,Energy2 s1,
				Energy m1,Energy m2,Energy m3);

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
   * @param F1 The form factor \f$F_1\f$.
   * @param F2 The form factor \f$F_2\f$.
   * @param F3 The form factor \f$F_3\f$.
   * @param F4 The form factor \f$F_4\f$.
   * @param F5 The form factor \f$F_5\f$.
   */
  virtual void calculateFormFactors(const int ichan,const int imode,
				    Energy2 q2,Energy2 s1,Energy2 s2,Energy2 s3,
				    Complex&F1,Complex&F2,Complex&F3,
				    Complex&F4,Complex&F5) const;

  /**
   * Calculate CLEO form factors for the current. Implements the form factors
   * described above.
   * @param imode The mode
   * @param ichan The phase space channel
   * @param q2 The scale \f$q^2\f$ for the current.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param F1 The form factor \f$F_1\f$.
   * @param F2 The form factor \f$F_2\f$.
   * @param F3 The form factor \f$F_3\f$.
   */
  void CLEOFormFactor(int imode,int ichan,Energy2 q2,Energy2 s1, Energy2 s2,
		      Energy2 s3,Complex & F1, Complex & F2, Complex & F3) const;

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
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  inline virtual void doinitrun();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ThreePionCLEOCurrent> initThreePionCLEOCurrent;

  /**
   * Private and non-existent assignment operator.
   */
  ThreePionCLEOCurrent & operator=(const ThreePionCLEOCurrent &);

private:

  /**
   * The \f$a_1\f$ running width
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @return The \f$a_1\f$ running width.
   */
  inline Energy a1width(Energy2 q2) const ;

  /**
   * Initialize the \f$a_1\f$ running width
   * @param iopt Initialization option (-1 full calculation, 0 set up the interpolation)
   */
  inline void inita1width(int iopt);
  
  /**
   * \f$a_1\f$ Breit-Wigner
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @return The Breit-Wigner
   */
  inline Complex a1BreitWigner(Energy2 q2) const;

  /**
   * The \f$\rho\f$ Breit-Wigner.
   * @param ires The \f$\rho\f$ multiplet.
   * @param q2 The scale, \f$q^2\f$.
   * @param icharge The charge of the \f$\rho\f$.
   * @return The Breit-Wigner
   */
  inline Complex rhoBreitWigner(int ires, Energy2 q2,int icharge) const;

  /**
   * Breit-Wigner for the \f$\sigma\f$.
   * @param q2 The scale, \f$q^2\f$.
   * @param icharge Whether the pions produced in the meson decay are
   * charged (0) or neutral (1) 
   * @return The Breit-Wigner
   */
  inline Complex sigmaBreitWigner(Energy2 q2,int icharge) const;
  
  /**
   * Breit-Wigner for the \f$f_0(1370)\f$.
   * @param q2 The scale, \f$q^2\f$.
   * @param icharge Whether the pions produced in the meson decay are
   * charged (0) or neutral (1) 
   * @return The Breit-Wigner
   */
  inline Complex f0BreitWigner(Energy2 q2,int icharge) const;

  /**
   * Breit-Wigner for the \f$f_2\f$.
   * @param q2 The scale, \f$q^2\f$.
   * @param icharge Whether the pions produced in the meson decay are
   * charged (0) or neutral (1) 
   * @return The Breit-Wigner
   */
  inline Complex f2BreitWigner(Energy2 q2,int icharge) const;

private:
  
  /**
   * Masses of the \f$\rho\f$ resonances.
   */
  vector<Energy> _rhomass;

  /**
   * Widths of the \f$\rho\f$ resonances.
   */
  vector<Energy> _rhowidth;

  /**
   * Momenta of the decay products for neutral \f$\rho\f$ decay.
   */
  vector<Energy> _prhocc;

  /**
   * Momenta of the decay products for charged \f$\rho\f$ decay.
   */
  vector<Energy> _prhoc0;

  /**
   * Mass of the \f$f_2\f$ resonance
   */
  Energy _f2mass;

  /**
   * Width of the \f$f_2\f$ resonance
   */
  Energy _f2width;

  /**
   * Momenta of the decay products for \f$f_2\f$ decay to charged pions.
   */
  Energy _pf2cc;

  /**
   * Momenta of the decay products for \f$f_2\f$ decay to neutral pions.
   */
  Energy _pf200;

  /**
   * Mass of the \f$f_0(1370)\f$ resonance
   */
  Energy _f0mass;

  /**
   * Width of the \f$f_0(1370)\f$ resonance
   */
  Energy _f0width;

  /**
   * Momenta of the decay products for \f$f_0(1370)\f$ decay to charged pions.
   */
  Energy _pf0cc;

  /**
   * Momenta of the decay products for \f$f_0(1370)\f$ decay to neutral pions.
   */
  Energy _pf000;

  /**
   * Mass of the \f$\sigma\f$ resonance
   */
  Energy _sigmamass;

  /**
   * Width of the \f$\sigma\f$ resonance
   */
  Energy _sigmawidth;

  /**
   * Momenta of the decay products for \f$\sigma\f$ decay to charged pions.
   */
  Energy _psigmacc;

  /**
   * Momenta of the decay products for \f$\sigma\f$ decay to neutral pions.
   */
  Energy _psigma00;

  /**
   * Mass of the neutral pion.
   */
  Energy _mpi0;

  /**
   * Mass of the charged pion.
   */
  Energy _mpic;

  /**
   * The \f$a_1\f$ mass
   */
  Energy _a1mass;

  /**
   * The \f$a_1\f$ width
   */
  Energy _a1width;

  /**
   * Mass of the \f$K^*\f$ resonace
   */
  Energy _mKstar;

  /**
   * Mass of the \f$K\f$ resonace
   */
  Energy _mK;

  /**
   * Coupling for the \f$KK^*\f$ term in the running width.
   */
  double _gammk;

  /**
   * pion decay constant
   */
  Energy _fpi;

  /**
   *  The prefactor
   */
  InvEnergy _fact;

  /**
   * Magnitude of the \f$p\f$-wave couplings of the rho resonance, \f$g^P_{\rho_k}\f$, 
   * (\f$\beta_{1,2}\f$ in the CLEO paper.)
   */
  vector<double> _rhomagP;

  /**
   * Phase of the \f$p\f$-wave couplings of the rho resonance, \f$g^P_{\rho_k}\f$,
   * (\f$\beta_{1,2}\f$ in the CLEO paper.)
   */
  vector<double> _rhophaseP;

  /**
   *\f$p\f$-wave couplings of the rho resonance, \f$g^P_{\rho_k}\f$, 
   * (\f$\beta_{1,2}\f$ in the CLEO paper.)
   */
  vector<Complex> _rhocoupP;

  /**
   * Magnitude of the \f$d\f$-wave couplings of the rho resonance, \f$g^D_{\rho_k}\f$,
   * (\f$\beta_{3,4}\f$ in the CLEO paper.)
   */
  vector<InvEnergy2> _rhomagD;

  /**
   * Phase of the \f$d\f$-wave couplings of the rho resonance, \f$g^D_{\rho_k}\f$, 
   * (\f$\beta_{3,4}\f$ in the CLEO paper.)
   */
  vector<double>_rhophaseD;

  /**
   * \f$d\f$-wave couplings of the rho resonance, \f$g^D_{\rho_k}\f$,
   * (\f$\beta_{3,4}\f$ in the CLEO paper.)
   */
  vector<complex<InvEnergy2> > _rhocoupD;

  /**
   * Magntiude of the coupling of the \f$f_2\f$ resonance, \f$g_{f_2}\f$,
   * (\f$\beta_5\f$ in the CLEO paper.)
   */
  InvEnergy2 _f2mag;

  /**
   * Phase of the coupling of the \f$f_2\f$ resonance, \f$g_{f_2}\f$,
   * (\f$\beta_5\f$ in the CLEO paper.)
   */
  double _f2phase;

  /**
   * Coupling of the \f$f_2\f$ resonance, \f$g_{f_2}\f$,
   * (\f$\beta_5\f$ in the CLEO paper.)
   */
  complex<InvEnergy2> _f2coup;

  /**
   * Magntiude of the coupling of the \f$f_0(1370)\f$ resonance, \f$g_{f_0}\f$,
   * (\f$\beta_6\f$ in the CLEO paper.)
   */
  double _f0mag;

  /**
   * Phase of the coupling of the \f$f_0(1370)\f$ resonance, \f$g_{f_0}\f$,
   * (\f$\beta_6\f$ in the CLEO paper.)
   */
  double _f0phase;

  /**
   * Coupling of the \f$f_0(1370)\f$ resonance, \f$g_{f_0}\f$,
   * (\f$\beta_6\f$ in the CLEO paper.)
   */
  Complex _f0coup;

  /**
   * Magntiude of the coupling of the \f$\sigma\f$ resonance, \f$g_\sigma\f$,
   * (\f$\beta_7\f$ in the CLEO paper.)
   */
  double _sigmamag;

  /**
   * Phase of the coupling of the \f$\sigma\f$ resonance, \f$g_\sigma\f$,
   * (\f$\beta_7\f$ in the CLEO paper.)
   */
  double _sigmaphase;

  /**
   * Coupling of the \f$\sigma\f$ resonance, \f$g_\sigma\f$,
   * (\f$\beta_7\f$ in the CLEO paper.)
   */
  Complex _sigmacoup;

  /**
   * use local values of the mass parameters
   */
  bool _localparameters;
  
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
  Interpolator *_a1runinter;

  /**
   * Initialize the running \f$a_1\f$ width.
   */
  bool _initializea1;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ThreePionCLEOCurrent.
 */
template <>
 struct BaseClassTrait<Herwig::ThreePionCLEOCurrent,1> {
  /** Typedef of the base class of ThreePionCLEOCurrent. */
  typedef Herwig::ThreeMesonCurrentBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ThreePionCLEOCurrent>
  : public ClassTraitsBase<Herwig::ThreePionCLEOCurrent> {
  /** Return the class name. */
  static string className() { return "Herwig++::ThreePionCLEOCurrent"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwWeakCurrents.so"; }

};

}
  
#include "CLHEP/GenericFunctions/AbsFunction.hh"
namespace Herwig {
using namespace Genfun;
using namespace ThePEG; 

/** \ingroup Decay
 *
 * Definitions of the functions to be integrated to give the running width
 * functions to return the matrix element for the a1 decay to be
 * integrated to give the a1 running width
 *
 * @see ThreeMesonDefaultCurrent
 */
class ThreePionCLEOa1MatrixElement: public Genfun::AbsFunction {

public:
  
  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;
  
  /**
   * Clone method
   */
   ThreePionCLEOa1MatrixElement *clone() const;

private:

  /**
   * Clone method
   */
  virtual AbsFunction *_clone() const;
    
public:
  
  /**
   * Constructor
   */
  ThreePionCLEOa1MatrixElement(int,ThreePionCLEOCurrentPtr);
  
  /**
   * Destructor
   */
  virtual ~ThreePionCLEOa1MatrixElement();

  /**
   * The  number of variables, in thsi case 7
   */  
  virtual unsigned int dimensionality() const ;     
  
  /**
   * Copy constructor
   */
  ThreePionCLEOa1MatrixElement(const ThreePionCLEOa1MatrixElement &right);
  
  /**
   * Retreive function value
   */
  virtual double operator ()(double) const {return 0.;}
  
  /**
   * Retreive function value
   */
  virtual double operator ()(const Argument & a) const ;
  
  /**
   * set the scale
   */
  inline void setQ2(Energy2);
  
  
private:
  
  /**
   * It is illegal to assign a function
   */
  const ThreePionCLEOa1MatrixElement & 
  operator=(const ThreePionCLEOa1MatrixElement &right);
  
  /**
   * the decayer
   * the integer for the mode
   */
private:
  Ptr<Herwig::ThreePionCLEOCurrent>::pointer _decayer;
  /**
   * the decayer
   * the integer for the mode
   */
  int _mode;
};

}

#include "ThreePionCLEOCurrent.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ThreePionCLEOCurrent.tcc"
#endif

#endif /* THEPEG_ThreePionCLEOCurrent_H */
