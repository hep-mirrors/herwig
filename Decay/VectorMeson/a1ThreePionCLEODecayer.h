// -*- C++ -*-
//
// a1ThreePionCLEODecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_a1ThreePionCLEODecayer_H
#define HERWIG_a1ThreePionCLEODecayer_H
//
// This is the declaration of the a1ThreePionCLEODecayer class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The  <code>a1ThreePionCLEODecayer</code> class is designed to implement the decay
 *  of the \f$a_1\f$ to three pions using the model of Phys.Rev.D61:012002,2000,
 *  (hep-ex/9902022) (CLEO) which was fitted to the one charged and two neutral pion
 *  channel for the charged \f$a_1\f$ decay in \f$\tau \to a_1 -> \pi\pi\pi\f$.
 *  The other modes are infered from this using isospin. This is a sophisticated model
 *  including the coupling of the \f$a_1\f$ to the \f$\rho\f$, \f$\rho(1450)\f$, 
 *  \f$f(1370)\f$ and \f$\sigma\f$ sigma mesons.
 *
 *  In this case the current is given by
 *  \f[\mathcal{M} = \epsilon_\mu
 *     \left[F_1(p_2-p_3)^\mu+F_2(p_3-p_1)^\mu+F_3(p_1-p_2)^\mu\right].\f]
 *
 *
 * The form factors for the \f$a_1^0 \to \pi^0 \pi^0 \pi^0\f$ mode are
 * \f[F_1=
 * \phantom{-}\frac23\left(g_\sigma B^S_\sigma(s_3)+g_{f_0}B^S_{f_0}(s_3)\right)
 *	   -\frac23\left(g_\sigma B^S_\sigma(s_2)+g_{f_0}B^S_{f_0}(s_2)\right)
 *         +g_{f_2}\left(\frac12(s_3-s_2)B^D_{f_2}(s_1)
 *   -\frac1{18}\frac{(4m_{\pi^0}^2-s_2)(q^2+s_2-m_{\pi^0}^2)}{s_2}B^D_{f_2}(s_2)
 *   +\frac1{18}\frac{(4m_{\pi^0}^2-s_3)(q^2-m_{\pi^0}^2+s_3)}{s_3}B^D_{f_2}(s_3)\right)
 *\f]
 *
 * \f[F_2=\phantom{-}\frac23(g_\sigma B^S_\sigma(s_3)+g_{f_0}B^S_{f_0}(s_3))
 *          -\frac23(g_\sigma B^S_\sigma(s_1)+g_{f_0}B^S_{f_0}(s_1))
 *          +g_{f_2}\left( \frac12(s_3-s1)B^D_{f_2}(s_2)
 *   -\frac1{18}\frac{(4m_{\pi^0}^2-s_1)(q^2+s_1-m_{\pi^0}^2)}{s_1}B^D_{f_2}(s_1)
 *   +\frac1{18}\frac{(4m_{\pi^0}^2-s_3)(q^2-m_{\pi^0}^2+s_3)}{s_3}B^D_{f_2}(s_3)\right)
 *\f]
 * \f[F_3=-\frac23(g_\sigma B^S_\sigma(s_1)+g_{f_0}B^S_{f_0}(s_1))
 *          +\frac23(g_\sigma B^S_\sigma(s_2)+g_{f_0}B^S_{f_0}(s_2))
 *          +g_{f_2}\left( \frac12(s_1-s_2)B^D_{f_2}(s_3)
 *  -\frac1{18}\frac{(4m_{\pi^0}^2-s_1)(q^2+s_1-m_{\pi^0}^2)}{s_1}B^D_{f_2}(s_1)
 *  +\frac1{18}\frac{(4m_{\pi^0}^2-s_2)(q^2+s_2-m_{\pi^0}^2)}{s_2}B^D_{f_2}(s_2)\right)
 *\f]
 *
 * The form factors for the \f$a_1^+ \to \pi^0 \pi^0 \pi^+\f$ mode are
 *
 * \f[F_1=\sum_k\left\{-\frac{g^P_{\rho_k}}3B_{\rho_k}^P(s_1)
 *          -g^D_{\rho_k}B_{\rho_k}^P(s_2)
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
 * 
 * The form factors for \f$a_1^0\to\pi^+\pi^-\pi^0\f$.
 *
 * \f[F_1=\sum_k\left\{g^P_{\rho_k}B_{\rho_k}^P(s_1)
 *    -\frac{g^D_{\rho_k}}3B_{\rho_k}^P(s_2)(s_3-m_{\pi^0}^2-s_1+m_{\pi^+}^2)\right\}
 * +\frac23\left(g_\sigma B^S_\sigma(s_3)+g_{f_0}B^S_{f_0}(s_3)\right)
 * +\frac{g_{f_2}}{18s_3}(q^2-m_{\pi^0}^2+s_3)(4m_{\pi^+}^2-s_3)B^D_{f_2}(s_3)\f]
 *
 * \f[F_2=\sum_k\left\{g^P_{\rho_k}B_{\rho_k}^P(s_2)
 *    -\frac{g^D_{\rho_k}}3B_{\rho_k}^P(s_1)(s_3-m_{\pi^0}^2-s_2+m_{\pi^+}^2)\right\}
 * +\frac23\left(g_\sigma B^S_\sigma(s_3)+g_{f_0}B^S_{f_0}(s_3)\right)
 * +\frac{g_{f_2}}{18s_3}(q^2-m_{\pi^0}^2+s_3)(4m_{\pi^+}^2-s_3)B^D_{f_2}(s_3)\f]
 *
 * \f[F_3=\sum_k
 *   g^D_{\rho_k}\left\{-\frac13B_{\rho_k}^P(s_1)(s_3-m_{\pi^0}^2-s_2+m_{\pi^+}^2)
 *                       +\frac13B_{\rho_k}^P(s_2)(s_3-m_{\pi^0}^2-s_1+m_{\pi^+}^2)
 *			\right\}
 * -\frac{g_{f_2}}2(s_1-s_2)B^D_{f_2}(s_3)\f]
 * 
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
    \f$B^L_Y(s_i) = \frac{m^2_Y}{m^2_Y-s_i)+im_Y\Gamma^{Y,L}(s_i)}\f$
 * where
 * \f$\Gamma^{Y,L}(s_i) = \Gamma^Y\left(\frac{p(s_i)}{p(M_Y}\right)^{2L+1}\frac{m_Y}{\sqrt{s_i}}\f$
 * \f$m_Y\f$ and \f$\Gamma^Y\f$ are the mass and width of the particle \f$Y\f$ 
 * respectively. \f$p(s_i)\f$ is the momentum of the outgoing pion in the 
 * rest frame of the resonanc \f$Y\f$.
 *
 * @see ThreePionCLEOCurrent
 * @see DecayIntegrator
 *
 */
class a1ThreePionCLEODecayer: public DecayIntegrator {
  
public:
  
  /**
   * Default constructor.
   */
  a1ThreePionCLEODecayer();

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
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const ParticleVector & decay, MEOption meopt) const;

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
  virtual IBPtr clone() const { return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this);}
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
  //@}
  
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<a1ThreePionCLEODecayer> inita1ThreePionCLEODecayer;
  
  /**
   * Private and non-existent assignment operator.
   */
  a1ThreePionCLEODecayer & operator=(const a1ThreePionCLEODecayer &);
  
private:

  /**
   * Breit wigner for the \f$\rho\f$, \f$B^P_{\rho_k}(q^2)\f$.
   * @param ires The \f$\rho\f$ multiplet to used.
   * @param q2 The scale, \f$q^2\f$.
   * @param icharge Which pion masses to use for the momentum calculation
   * @return The Breit-Wigner
   */
  Complex rhoBreitWigner(int ires, Energy2 q2,int icharge) const {
    Energy q=sqrt(q2);
    Complex ii(0.,1.);
    double ratio = icharge==0 ? 
      Kinematics::pstarTwoBodyDecay(q,_mpic,_mpic)/_prhocc[ires] :
      Kinematics::pstarTwoBodyDecay(q,_mpic,_mpi0)/_prhoc0[ires];
    Energy gamrun=_rhowidth[ires]*pow(ratio,3)*_rhomass[ires]/q;
    return sqr(_rhomass[ires])
      /(sqr(_rhomass[ires])-q2-ii*_rhomass[ires]*gamrun);
  }

  /**
   * Breit wigner for the \f$\sigma\f$, \f$B^S_\sigma(q^2)\f$.
   * @param q2 The scale, \f$q^2\f$.
   * @param icharge Which pion masses to use for the momentum calculation
   * @return The Breit-Wigner
   */
  Complex sigmaBreitWigner(Energy2 q2,int icharge) const {
    Energy q=sqrt(q2);
    Complex ii(0.,1.);
    double ratio = icharge==0 ? 
      Kinematics::pstarTwoBodyDecay(q,_mpic,_mpic)/_psigmacc :
      Kinematics::pstarTwoBodyDecay(q,_mpi0,_mpi0)/_psigma00;
    Energy gamrun=_sigmawidth*ratio*_sigmamass/q;
    return sqr(_sigmamass)/(sqr(_sigmamass)-q2-ii*_sigmamass*gamrun);
  }
  
  /**
   * Breit wigner for the \f$f_0(1370)\f$, \f$B^S_{f_0}(q^2)\f$.
   * @param q2 The scale, \f$q^2\f$.
   * @param icharge Which pion masses to use for the momentum calculation
   * @return The Breit-Wigner
   */
  Complex f0BreitWigner(Energy2 q2,int icharge) const {
    Energy q=sqrt(q2);
    Complex ii(0.,1.);
    double ratio = icharge==0 ? 
      Kinematics::pstarTwoBodyDecay(q,_mpic,_mpic)/_pf0cc :
      Kinematics::pstarTwoBodyDecay(q,_mpi0,_mpi0)/_pf000;
    Energy gamrun=_f0width*ratio*_f0mass/q;
    return sqr(_f0mass)/(sqr(_f0mass)-q2-ii*_f0mass*gamrun);
  }
  
  /**
   * Breit wigner for the \f$f_2\f$, \f$B^D_{f_2}(q^2)\f$.
   * @param q2 The scale, \f$q^2\f$.
   * @param icharge Which pion masses to use for the momentum calculation
   * @return The Breit-Wigner
   */
  Complex f2BreitWigner(Energy2 q2,int icharge) const {
    Energy q=sqrt(q2);
    Complex ii(0.,1.);
    double ratio = icharge==0 ?
      Kinematics::pstarTwoBodyDecay(q,_mpic,_mpic)/_pf2cc :
      Kinematics::pstarTwoBodyDecay(q,_mpi0,_mpi0)/_pf200;
    Energy gamrun=_f2width*pow(ratio,5)*_f2mass/q;
    return sqr(_f2mass)/(sqr(_f2mass)-q2-ii*_f2mass*gamrun);
  }
  
  /**
   * Calculate the form factors
   * @param iopt The mode being calculated in the order given above
   * @param ichan The phase space channel in the order given in the doinit member.
   * @param q2 The sacale \f$q^2\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param F1 The form factor \f$F_1\f$.
   * @param F2 The form factor \f$F_2\f$.
   * @param F3 The form factor \f$F_3\f$.
   * 
   */
  void formFactors(int iopt,int ichan,Energy2 q2,Energy2 s1,Energy2 s2,
		   Energy2 s3,
		   complex<InvEnergy> & F1,
		   complex<InvEnergy> & F2,
		   complex<InvEnergy> & F3) const;

private:

  /**
   * Masses of the rho resonaces
   */
  vector<Energy> _rhomass;

  /**
   * Widths of the rho resonaces
   */
  vector<Energy> _rhowidth;

  /**
   * Momentum of the particles produced in charged rho decay
   */
  vector<Energy> _prhocc;

  /**
   * Momentum of the particles produced in neutral rho decay
   */
  vector<Energy> _prhoc0;

  /**
   * Mass of the \f$f_2\f$.
   */
  Energy _f2mass;

  /**
   * Width of the \f$f_2\f$.
   */
  Energy _f2width;

  /**
   * Momentum for the decay of the \f$f_2\f$ to two charged pions.
   */
  Energy _pf2cc;

  /**
   * Momentum for the decay of the \f$f_2\f$ to two neutral pions.
   */
  Energy _pf200;

  /**
   * Mass of the \f$f_0(1370)\f$.
   */
  Energy _f0mass;

  /**
   * Width of the \f$f_0(1370)\f$.
   */
  Energy _f0width;

  /**
   * Momentum for the decay of the \f$f_0(1370)\f$ to two charged pions.
   */
  Energy _pf0cc;

  /**
   * Momentum for the decay of the \f$f_0(1370)\f$ to two neutral pions.
   */
  Energy _pf000;

  /**
   * Mass of the \f$\sigma\f$ meson.
   */
  Energy _sigmamass;

  /**
   * Width of the \f$\sigma\f$ meson.
   */
  Energy _sigmawidth;

  /**
   * Momentum for the decay of the \f$\sigma\f$ to two charged pions.
   */
  Energy _psigmacc;

  /**
   * Momentum for the decay of the \f$\sigma\f$ to two neutral pions.
   */
  Energy _psigma00;

  /**
   * Mass of the neutral pion
   */
  Energy _mpi0;

  /**
   * Mass of the charged pion
   */
  Energy _mpic;

  /**
   * overall coupling for the decay
   */
  InvEnergy _coupling;

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
   * Use local values of the mass parameters
   */
  bool _localparameters;
  
  /**
   * Weights for the channels for the zero charged pion channel.
   */
  mutable vector<double> _zerowgts;
  
  /**
   * Weights for the channels for the one charged pion channel.
   */
  mutable vector<double> _onewgts;
  
  /**
   * Weights for the channels for the two charged pion channel.
   */
  mutable vector<double> _twowgts;
  
  /**
   * Weights for the channels for the three charged pion channel.
   */
  mutable vector<double> _threewgts;

  /**
   * Maximum weight for the zero charged pion channel.
   */
  mutable double _zeromax;

  /**
   * Maximum weight for the one charged pion channel.
   */
  mutable double _onemax;

  /**
   * Maximum weight for the two charged pion channel.
   */
  mutable double _twomax;

  /**
   * Maximum weight for the three charged pion channel.
   */
  mutable double _threemax;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization vectors
   */
  mutable vector<Helicity::LorentzPolarizationVector> _vectors;
};
  
}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */
  
/**
 * The following template specialization informs ThePEG about the
 * base class of a1ThreePionCLEODecayer.
 */
template <>
struct BaseClassTrait<Herwig::a1ThreePionCLEODecayer,1> {
  /** Typedef of the base class of a1ThreePionCLEODecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};
  
template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
struct ClassTraits<Herwig::a1ThreePionCLEODecayer>
  : public ClassTraitsBase<Herwig::a1ThreePionCLEODecayer> {
  /** Return the class name. */
  static string className() { return "Herwig::a1ThreePionCLEODecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwVMDecay.so"; }
  
};

/** @endcond */
  
}

#endif /* HERWIG_a1ThreePionCLEODecayer_H */
