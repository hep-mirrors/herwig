// -*- C++ -*-
//
// a1ThreePionDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_a1ThreePionDecayer_H
#define HERWIG_a1ThreePionDecayer_H
//
// This is the declaration of the a1ThreePionDecayer class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 *  The  <code>a1ThreePionDecayer</code> class is designed to implement the decay
 *  of the a_1 to three pions. The model used is one based on the Novosibirsk
 *  four pion current used in TAUOLA. 
 *
 *  - The matrix element for the decay \f$a_1^0\to\pi^0\pi^0\pi^0\f$ is
 * 
 *    \f[J^\mu =\frac{F_{a_1}(q^2)g}{qM^2_{\rho_0}}\epsilon_\mu\left[
 * q^2z \sum_i p^\mu_i \frac1{D_\sigma(s_i)}\right]\f]
 *  
 *  - The matrix element for the decay \f$a_1^+\to\pi^0\pi^0\pi^+\f$ is
 *
 *    \f[J^\mu =\frac{F_{a_1}(q^2)g}{qM^2_{\rho_0}}\epsilon_\mu\left[
 *  q^2z\frac1{D_\sigma(s_3)} p_3^\mu 
 *  +\sum_k g_{\rho_k}\left\{
 *      \frac1{D_{\rho_k}(s_1)}\left(p_0\cdot p_3p_2^\mu-p_0\cdot p_2p_3^\mu\right)
 *     +\frac1{D_{\rho_k}(s_2)}\left(p_0\cdot p_3p_1^\mu-p_0\cdot p_1p_3^\mu\right) \right\}
 * \right]\f]
 *
 *  - The matrix element for the decay \f$a_10\to\pi^+\pi^-\pi^0\f$ is
 *
 *   \f[J^\mu =\frac{F_{a_1}(q^2)g}{qM^2_{\rho_0}}\left[ 
 * q^2z\frac1{D_\sigma(s_3)}p_3^\mu
 *  +\sum_k g_{\rho_k}\left\{
 *      \frac1{D_{\rho_k}(s_1)}\left(p_0\cdot p_3p_2^\mu-p_0\cdot p_2p_3^\mu\right)
 *     +\frac1{D_{\rho_k}(s_2)}\left(p_0\cdot p_3p_1^\mu-p_0\cdot p_1p_3^\mu\right)\right\}
 * \right]\f]
 *
 *  - The current for the decay \f$a_1^+\to\pi^+\pi^+\pi^-\f$ is
 *
 *   \f[J^\mu = \frac{F_{a_1}(q^2)g}{qM^2_{\rho_0}}\left[
 * q^2z\left(\frac1{D_\sigma(s_1)}p_1^\mu+\frac1{D_\sigma(s_2)}p_2^\mu\right)
 *  -\sum_k g_{\rho_k}\left\{
 *      \frac1{D_{\rho_k}(s_1)}\left(p_0\cdot p_3p_2^\mu-p_0\cdot p_2p_3^\mu\right)
 *     +\frac1{D_{\rho_k}(s_2)}\left(p_0\cdot p_3p_1^\mu-p_0\cdot p_1p_3^\mu\right)\right\}
 * \right]\f]
 *
 *  *  The denominator factor is
 *  \f[D_\sigma(q^2) = q^2-M^2+iM\Gamma\frac{g_\sigma(q^2)}{g_\sigma(M^2)}\f] 
 *  where \f$g(s) = \left(1-4\frac{m_{\pi}^2}{s}\right)\f$ for the \f$\sigma\f$ meson,
 *  and
 * \f[D_{\rho_k}(q^2) = q^2-M^2_{\rho_k}-M_{\rho_k}\Gamma_{\rho_k}dm(q^2)
 *   +iM_{\rho_k}\Gamma_{\rho_k}\frac{g_{\rho_k}(q^2)}{g_{\rho_k}(M^2)}
 * \f]
 *  for the \f$\rho\f$.
 *  
 *  The propagator factors are normalized such that \f$D(0)=-1\f$.
 *
 *  Here
 *  \f[dm(q^2) = \frac1{g_{\rho_k}(M^2_{\rho_k}}\left(h_{\rho_k}(q^2)-h_{\rho_k}(M^2_{\rho_k}
 *   -\left.(q^2-M^2_{\rho_k})\frac{dh_{\rho_k}(q^2)}{dq^2}\right|_{q^2=M^2_{\rho_k}}
 *   \right)\f]
 *  where
 * \f[h_{\rho_k}(q^2) = 
 *   \left\{ \begin{array}{cc}
 *   \frac{\sqrt{1-\frac{4m_\pi^2}{q^2}}\ln\left(\frac{1+\sqrt{1-\frac{4m_\pi^2}{q^2}}}{1-\sqrt{1-\frac{4m_\pi^2}{q^2}}}\right)(q^2-4m_\pi^2)}\pi  & {\rm for\ } q^2>4m_\pi^2 \\
 *   -8\frac{m_\pi^2}{\pi} & {\rm for\ } q^2=0\,{\rm GeV^2} \\
 *   0                     & {\rm otherwise} 
 *   \end{array}\right.\f]  
 *
 *  The \f$a_1\f$ form factor for the off-shell \f$a_1\f$ is given by
 * \f[F_{a_1}(q^2) = \frac{\left(1+m^2_{a_1}/\Lambda^2\right)}
 *                        {\left(1+      q^2/\Lambda^2\right)}.\f]
 *   
 *  The masses and couplings are
 * - \f$z\f$ is the relative coupling of the \f$\sigma\f$.
 * - \f$m_\pi\f$ is the mass of the pion.
 * - \f$g_{\rho_k}\f$ is the coupling of the \f$k\f$th \f$\rho\f$ multiplet.
 * - \f$M_{\rho_k}\f$ the mass of the \f$k\f$th \f$\rho\f$ resonance.
 * - \f$\Gamma_{\rho_k}\f$ the width of the \f$k\f$th \f$\rho\f$ resonance.
 * - \f$M_\sigma\f$ the mass of the \f$\sigma\f$ meson.
 * - \f$\Gamma_\sigma\f$ the width of the \f$\sigma\f$ meson.
 * - \f$m_{a_1}\f$ the mass of the \f$a_1\f$ meson.
 * - \f$\Lambda^2\f$ the mass parameter for the \f$a_1\f$ form factor.
 * @see FourPionNovosibirskCurrent
 * @see DecayIntegrator
 * 
 */
class a1ThreePionDecayer: public DecayIntegrator {
  
public:
  
  /**
   * Default constructor.
   */
  a1ThreePionDecayer();

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
  //@}

private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<a1ThreePionDecayer> inita1ThreePionDecayer;
  
  /**
   * Private and non-existent assignment operator.
   */
  a1ThreePionDecayer & operator=(const a1ThreePionDecayer &) = delete;
  
private:
  
  /**
   * Breit-wigner for the \f$\sigma\f$, this is \f$\frac1{D_\sigma(q^2)}\f$.
   * @param q2 The scale, \f$q^2\f$.
   * @return The Breit-Wigner
   */
  Complex sigmaBreitWigner(Energy2 q2) const {
    Energy q=sqrt(q2);
    Energy width=_sigmawidth*Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi)/_psigma;
    Energy2 msigma2=_sigmamass*_sigmamass;
    Complex ii(0.,1.);
    complex<Energy2> denom = q>2.*_mpi ? q2-msigma2+ii*msigma2*width/q :
      q2-msigma2;
    return msigma2/denom;
  }
  
  /**
   * The \f$a_1\f$ form factor, \f$F_{a_1}(q^2)\f$
   * @param q2 The scale, \f$q^2\f$.
   * @return The form factor.
   */
  double a1FormFactor(Energy2 q2) const {
    return (1.+_a1mass2/_lambda2)/(1.+q2/_lambda2);
  }

  /**
   * Breit-Wigner for the \f$\rho\f$, this is  \f$\frac1{D_{\rho_k}(q^2)}\f$.
   * @param q2 The scale, \f$q^2\f$.
   * @param ires The \f$\rho\f$ multiplet.
   * @return The Breit-Wigner
   */
  Complex rhoBreitWigner(Energy2 q2,int ires) const {
    Energy q=sqrt(q2);
    Energy2 grhom = 8.*_prho[ires]*_prho[ires]*_prho[ires]/_rhomass[ires];
    complex<Energy2> denom;
    Complex ii(0.,1.);
    if(q2<4.*_mpi2) {
      denom=q2-_rhomass[ires]*_rhomass[ires]-_rhowidth[ires]*_rhomass[ires]*
	(hFunction(q)-_hm2[ires]-(q2-_rhomass[ires]*_rhomass[ires])*_dhdq2m2[ires])
	/grhom;
    }
    else {
      Energy pcm=2.*Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
      Energy2 grho = pcm*pcm*pcm/q;
      denom=q2-_rhomass[ires]*_rhomass[ires]
	-_rhowidth[ires]*_rhomass[ires]*
	(hFunction(q)-_hm2[ires]-(q2-_rhomass[ires]*_rhomass[ires])*_dhdq2m2[ires])/grhom
	+ii*_rhomass[ires]*_rhowidth[ires]*grho/grhom;
    }
    return _rhoD[ires]/denom;
  }

  /**
   *  Normalisation factor for the \f$\rho\f$ propagator to ensure \f$D(0)=-1\f$.
   * @param ires The \f$\rho\f$ multiplet.
   * @return The normalisation factor.
   */
  Energy2 DParameter(int ires) const {
    Energy2 grhom = 8.*_prho[ires]*_prho[ires]*_prho[ires]/_rhomass[ires];
    return _rhomass[ires]*_rhomass[ires]+_rhowidth[ires]*_rhomass[ires]*
      (hFunction(ZERO)-_hm2[ires]+sqr(_rhomass[ires])*_dhdq2m2[ires])/grhom;
  }

  /**
   * The \f$\frac{dh}{dq^2}\f$ function in the rho propagator evaluated at \f$q^2=m^2\f$.
   * @param ires The \f$\rho\f$ resonance for the function
   * @return \f$\frac{dh}{dq^2}\f$ evaluated at \f$q^2=m^2\f$.
   */
  double dhdq2Parameter(int ires) const  {
    Energy2 mrho2(sqr(_rhomass[ires]));
    double root = sqrt(1.-4.*_mpi2/mrho2);
    using Constants::pi;
    return root/pi*(root+(1.+2*_mpi2/mrho2)*log((1+root)/(1-root)));
  }
  
  /**
   * The \f$h(q^2)\f$ function in the \f$\rho\f$ propagator.
   * @param q The scale, \f$q\f$.
   * @return The function \f$h(q^2)\f$.
   */
  Energy2 hFunction(const Energy q) const  {
    static const Energy2 eps(0.01*MeV2);
    Energy2 q2=sqr(q), output;
    double root = sqrt(1.-4.*_mpi2/q2);
    if(q2>4*_mpi2) {
      output=root*log((1.+root)/(1.-root))*(q2-4*_mpi2)/Constants::pi;
    }
    else if(q2>eps) output=ZERO;
    else            output=-8.*_mpi2/Constants::pi;
    return output;
  }

  /**
   *  Momentum Function
   */
  Energy4 lambda(Energy2 a, Energy2 b, Energy2 c) const {
    return sqr(a)+sqr(b)+sqr(c)-2.*a*b-2.*a*c-2.*b*c;
  }
  
private:

  /**
   * Mass of the rho resonances
   */
  vector<Energy> _rhomass;

  /**
   * Width of the rho resonaces
   */
  vector<Energy> _rhowidth;

  /**
   * Momentum of the pions produced in the \f$\rho\f$ decay.
   */
  vector<Energy> _prho;

  /**
   * The function \f$h(q^2)\f$ evaluated at \f$q^2=M^2_{\rho_k}\f$
   */
  vector<Energy2> _hm2;

  /**
   * The normalization factor for the \f$\rho_k\f$ propagator factor.
   */
  vector<Energy2> _rhoD;

  /**
   * The \f$\frac{dh}{dq^2}\f$ function in the rho propagator evaluated at \f$q^2=m^2\f$
   * for the different \f$\rho\f$ multiplets.
   */
  vector<double> _dhdq2m2;

  /**
   * The mass of the \f$\sigma\f$ meson. 
   */
  Energy _sigmamass;

  /**
   * The width of the \f$\sigma\f$ meson. 
   */
  Energy _sigmawidth;

  /**
   * The momenta of the pions produced in the \f$\sigma\f$ meson decay. 
   */
  Energy _psigma;

  /**
   * The mass of the pion, \f$m_\pi\f$.
   */
  Energy _mpi;

  /**
   * The mass of the pion, \f$m^2_\pi\f$.
   */
  Energy2 _mpi2;

  /**
   * The \f$\Lambda^2\f$ parameter for the \f$a_1\f$ form factor.
   */
  Energy2 _lambda2;
  
  /**
   * The mass squared of the \f$a_1\f$ meson, \f$m_{a_1}^2\f$.
   */
  Energy2 _a1mass2;

  /**
   * The \f$z\f$ coupling for the \f$\sigma\f$ resonance.
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
   * \f$g_{\rho_k}\f$ is the coupling of the \f$k\f$ th \f$\rho\f$ multiplet.
   */
  vector<Complex> _rhocoupling;

  /**
   *  Magnitude of the rho coupling
   */
  vector<double> _rhomag;

  /**
   *  Phase of the rho coupling
   */
  vector<double> _rhophase;

  /**
   * The overall coupling for the decay.
   */
  double _coupling;

  /**
   * use local values of the mass parameters
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
 * base class of a1ThreePionDecayer.
 */ 
template <>
struct BaseClassTrait<Herwig::a1ThreePionDecayer,1> {
  /** Typedef of the base class of a1ThreePionDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::a1ThreePionDecayer>
  : public ClassTraitsBase<Herwig::a1ThreePionDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig::a1ThreePionDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwVMDecay.so"; }
  
};

/** @endcond */
  
}

#endif /* HERWIG_a1ThreePionDecayer_H */
