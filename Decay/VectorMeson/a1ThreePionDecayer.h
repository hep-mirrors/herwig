// -*- C++ -*-
#ifndef HERWIG_a1ThreePionDecayer_H
#define HERWIG_a1ThreePionDecayer_H
//
// This is the declaration of the a1ThreePionDecayer class.
//
#include "VectorMesonDecayerBase.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig++/Utilities/Kinematics.h"
// #include "a1ThreePionDecayer.fh"
// #include "a1ThreePionDecayer.xh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 *  The  <code>a1ThreePionDecayer</code> class is designed to implement the decay
 *  of the a_1 to three pions. The model used is one based on the Novosibirsk
 *  four pion current used in TAUOLA. 
 *
 *  - The current for the decay \f$a_1^0\to\pi^0\pi^0\pi^0\f$ is
 * 
 *    \f[J^\mu =\frac{F_{a_1}(q^2)g}{qM^2_{\rho_0}}\left[
 * q^2z \sum_i p^\mu_i \frac1{D_\sigma(s_i)}\right]\f]
 *  
 *  - The current for the decay \f$a_1^+\to\pi^0\pi^0\pi^+\f$ is
 *
 *    \f[J^\mu =\frac{F_{a_1}(q^2)g}{qM^2_{\rho_0}}\left[
 *  q^2z\frac1{D_\sigma(s_3)} p_3^\mu 
 *  +\sum_k g_{\rho_k}\left\{
 *      \frac1{D_{\rho_k}(s_1)}\left(p_0\cdot p_3p_2^\mu-p_0\cdot p_2p_3^\mu\right)
 *     +\frac1{D_{\rho_k}(s_2)}\left(p_0\cdot p_3p_1^\mu-p_0\cdot p_1p_3^\mu\right) \right\}
 * \right]\f]
 *
 *  - The current for the decay \f$a_10\to\pi^+\pi^-\pi^0\f$ is
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
 * @see VectorMesonDecayerBase
 * 
 */
class a1ThreePionDecayer: public VectorMesonDecayerBase {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline a1ThreePionDecayer();

  /**
   * Copy-constructor.
   */
  inline a1ThreePionDecayer(const a1ThreePionDecayer &);

  /**
   * Standard ctors and dtor.
   */
  virtual ~a1ThreePionDecayer();
  //@}
  
public:
  
  /**
   * Accept member which is called at initialization to see if this Decayer can
   * handle a given decay mode. This version checks the particles against the 
   * list of allowed incoming  and outgoing mesons.
   * @param dm The DecayMode
   * @return Whether the mode can be handled.
   */
  virtual bool accept(const DecayMode & dm) const;
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. This version uses PDG codes to
   * work which mode is being simulated and the generate member of the 
   * DecayIntegrator class for the phase-space  generation.
   * @param dm The DecayMode
   * @param part The Particle instant being decayed.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & part) const;
  
  /**
   * The hadronic current. This returns the current 
   * described above.
   * @param vertex Construct the information for spin correlations.
   * @param ichan The phase-space channel to calculate the current for.
   * @param inpart The decaying particle
   * @param outpart The decay products
   * @return The hadronic currents for the decay.
   */
  virtual vector<LorentzPolarizationVector> 
  decayCurrent(const bool vertex, const int ichan,const Particle & inpart, 
	       const ParticleVector & outpart) const;
  
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
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<a1ThreePionDecayer> inita1ThreePionDecayer;
  
  /**
   * Private and non-existent assignment operator.
   */
  a1ThreePionDecayer & operator=(const a1ThreePionDecayer &);
  
private:
  
  /**
   * Breit-wigner for the \f$\sigma\f$, this is \f$\frac1{D_\sigma(q^2)}\f$.
   * @param q2 The scale, \f$q^2\f$.
   * @return The Breit-Wigner
   */
  inline complex<InvEnergy2> sigmaBreitWigner(Energy2 q2) const;
  
  /**
   * The \f$a_1\f$ form factor, \f$F_{a_1}(q^2)\f$
   * @param q2 The scale, \f$q^2\f$.
   * @return The form factor.
   */
  inline double a1FormFactor(Energy2 q2) const;

  /**
   * Breit-Wigner for the \f$\rho\f$, this is  \f$\frac1{D_{\rho_k}(q^2)}\f$.
   * @param q2 The scale, \f$q^2\f$.
   * @param ires The \f$\rho\f$ multiplet.
   * @return The Breit-Wigner
   */
  inline Complex rhoBreitWigner(Energy2 q2,int ires) const;

  /**
   *  Normalisation factor for the \f$\rho\f$ propagator to ensure \f$D(0)=-1\f$.
   * @param ires The \f$\rho\f$ multiplet.
   * @return The normalisation factor.
   */
  inline Energy2 DParameter(int ires) const;

  /**
   * The \f$\frac{dh}{dq^2}\f$ function in the rho propagator evaluated at \f$q^2=m^2\f$.
   * @param ires The \f$\rho\f$ resonance for the function
   * @return \f$\frac{dh}{dq^2}\f$ evaluated at \f$q^2=m^2\f$.
   */
  inline double dhdq2Parameter(int ires) const ;
  
  /**
   * The \f$h(q^2)\f$ function in the \f$\rho\f$ propagator.
   * @param q2 The scale, \f$q^2\f$.
   * @return The function \f$h(q^2)\f$.
   */
  inline Energy2 hFunction(const Energy q2) const ;
  
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
   * The mass of the \f$a_1\f$ meson, \f$m_{a_1}\f$.
   */
  Energy _a1mass2;

  /**
   * The \f$z\f$ coupling for the \f$\sigma\f$ resonance.
   */
  Complex _zsigma;

  /**
   * \f$g_{\rho_k}\f$ is the coupling of the \f$k\f$th \f$\rho\f$ multiplet.
   */
  vector<Complex> _rhocoupling;

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

};
  
}


namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of a1ThreePionDecayer.
   */
  template <>
  struct BaseClassTrait<Herwig::a1ThreePionDecayer,1> {
    /** Typedef of the base class of a1ThreePionDecayer. */
    typedef Herwig::VectorMesonDecayerBase NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::a1ThreePionDecayer>
    : public ClassTraitsBase<Herwig::a1ThreePionDecayer> {
    /** Return the class name.*/
    static string className() { return "Herwig++::a1ThreePionDecayer"; }
    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwVMDecay.so"; }

  };
  
}

#include "a1ThreePionDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "a1ThreePionDecayer.tcc"
#endif

#endif /* HERWIG_a1ThreePionDecayer_H */
