// -*- C++ -*-
#ifndef HERWIG_VectorMeson3PionDecayer_H
#define HERWIG_VectorMeson3PionDecayer_H
// This is the declaration of the VectorMeson3PionDecayer class.

#include "VectorMesonDecayerBase.h"
#include "Herwig++/Utilities/Kinematics.h"
// #include "VectorMeson3PionDecayer.fh"
// #include "VectorMeson3PionDecayer.xh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>VectorMeson3PionDecayer</code> class is designed to perform the
 *  decay of an \f$I=0\f$ meson to three pions via rho mesons including the
 *  option of higher rho resonaces and a constant term. It is mainly intended for
 *  the decays
 *
 *   - \f$\omega \to \pi^+\pi^-\pi^0\f$
 *   - \f$\phi   \to \pi^+\pi^-\pi^0\f$
 *
 *  The default for the \f$omega\f$ is to only include the
 *  contributions of the \f$\rho(770)\f$
 *  without a constant term, whereas for the \f$\phi\f$ the parameters
 *  from hep-ex/0303016 (KLOE) which includes the \f$\rho(770)\f$ 
 *  and a constant term to represent the effects
 *  of the higher rho resonances. (The KLOE paper also included a omega contribution
 *  but this is assumed to be non-resonant.)
 *
 *  The form of the current is taken to be
 *
 *  \f[ J^\mu = g\epsilon^{\mu\alpha\alpha\beta\gamma}
 *                      p_{+\alpha}p_{-\beta}p_{0\gamma}
 *      \left[a_de^{i\phi_d} +\sum_k a_{\rho_k}e^{i\phi_{\rho_k}}
 *   \left\{
 *  \frac{M^2_{\rho_k^+}}{m^2_{0+}-M^2_{\rho_k^+}
 *           +im_{0+}\Gamma_{\rho_k^+}(m^2_{0+})}
 *  +\frac{M^2_{\rho_k^-}}{m^2_{0-}-M^2_{\rho_k^-}
 *           +im_{0-}\Gamma_{\rho_k^-}(m^2_{0-})}
 *  +\frac{M^2_{\rho_k^-}}{m^2_{+-}-M^2_{\rho_k^+}
 *           +im_{+-}\Gamma_{\rho_k^0}(m^2_{+-})}\right\}
 *  \right],\f]
 *
 *  where \f$p_{+,-,0}\f$ are the momenta of the positively charged, negatively charged
 *  and neutral pions respectively, \f$m^2_{ij}=(p_i+p_j)^2\f$, \f$M_{\rho_k}\f$ is
 *  the mass of the \f$k\f$th \f$\rho\f$ resonance, \f$g\f$ is the overall coupling and 
 * \f[\Gamma{\rho_k}(m^2) = \Gamma_{\rho_k}\left(\frac{p_\pi(m^2)}{p_\pi(M_{\rho_k}^2)}\right)^2
 *    \left(\frac{M_{\rho_k}^2}{m^2}\right) \f]
 *  where \f$p_\pi\f$ is the pion momentum in the \f$\rho\f$ rest frame
 *  and \f$\Gamma_{\rho_k}\f$ is the width. In practice the couplings of the different
 *  terms are measured relative to the coupling of the lightest \f$\rho\f$ multiplet
 *   which is assumed to be one.
 *
 *
 *  To allow the easy addition of further modes the parameters for additional modes
 *  can be set. The following must be specified for each mode
 *
 *   - Incoming       - the PDG code for the incoming particle
 *   - Coupling       - the overall coupling for the decay, \f$g\f$.
 *   - DirectCoupling - the relative coupling for the direct term, \f$a_d\f$.
 *   - DirectPhase    - the phase of the coupling for the direct term, \f$\phi_d\f$.
 *   - Rho2Coupling   - the relative coupling for the second rho multiplet, 
 *     \f$a_{\rho_2}\f$.
 *   - Rho2Phase      - the phase of the coupling for the second rho multiplet, 
 *     \f$\phi_{\rho_2}\f$.
 *   - Rho3Coupling   - the relative coupling for the third rho multiplet, 
 *     \f$a_{\rho_3}\f$.
 *   - Rho3Phase      - the phase of the coupling for the third rho multiplet, 
 *     \f$\phi_{\rho_3}\f$.
 *   - MaximumWeight  - the maximum weight for the integration of the channel
 *   - Rho1Weight     - the weight for the first rho multiplet in the 
 *                      multichannel integration
 *   - Rho2Weight     - the weight for the second rho multiplet in the
 *                      multichannel integration
 *   - Rho3Weight     - the weight for the third rho multiplet in the
 *                      multichannel integration
 *   - Rho1Mass       - mass  of the first  rho multiplet, \f$M_{\rho_1}\f$.
 *   - Rho2Mass       - mass  of the second rho multiplet, \f$M_{\rho_2}\f$.
 *   - Rho3Mass       - mass  of the third  rho multiplet, \f$M_{\rho_3}\f$.
 *   - Rho1Width      - width of the first  rho multiplet, \f$\Gamma_{\rho_1}\f$.
 *   - Rho2Width      - width of the second rho multiplet, \f$\Gamma_{\rho_2}\f$.
 *   - Rho3Width      - width of the third  rho multiplet, \f$\Gamma_{\rho_3}\f$.
 *
 *
 * @see VectorMesonDecayerBase
 * 
 *  \author Peter Richardson
 */
class VectorMeson3PionDecayer: public VectorMesonDecayerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline VectorMeson3PionDecayer();

  /**
   * Copy-constructor.
   */
  inline VectorMeson3PionDecayer(const VectorMeson3PionDecayer &);

  /**
   * Destructor.
   */
  virtual ~VectorMeson3PionDecayer();
  //@}

public:


  /**
   * The hadronic current. This returns the current 
   *  described above.
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

  /**
   * Accept member which is called at initialization to see if this Decayer can
   * handle a given decay mode. This version checks that the outgoing pions are correct
   * and the incoming vector meson is one of those allowed.
   * @param dm The DecayMode
   * @return Whether the mode can be handled.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. This version uses PDG codes to
   * work which incoming meson is being simulated and the generate member of the 
   * DecayIntegrator class for the phase-space  generation.
   * @param dm The DecayMode
   * @param part The Particle instant being decayed.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & part) const;

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
  virtual double threeBodyMatrixElement(int imode,Energy2 q2, Energy2 s3,Energy2 s2,
					Energy2 s1,Energy m1,Energy m2,Energy m3);
  

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
  static ClassDescription<VectorMeson3PionDecayer> initVectorMeson3PionDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  VectorMeson3PionDecayer & operator=(const VectorMeson3PionDecayer &);

private:

  /**
   * vector storing the decaying particles for the different modes
   */
  vector<double> _incoming;

  /**
   * the overall coupling for the decay, \f$g\f$.
   */
  vector<double> _coupling;

  /**
   * relative coupling for the direct term, \f$a_d\f$.
   */
  vector<double> _directcoupling;

  /**
   * relative phase for the direct term, \f$\phi_d\f$.
   */
  vector<double> _directphase;

  /**
   * relative coupling for the second rho multiplet, \f$a_{\rho_2}\f$.
   */
  vector<double> _rho2coupling;

  /**
   * relative phase for the second rho multiplet, \f$\phi_{\rho_2}\f$.
   */
  vector<double> _rho2phase;

  /**
   * relative coupling for the second rho multiplet, \f$a_{\rho_3}\f$.
   */
  vector<double> _rho3coupling;

  /**
   * relative phase for the second rho multiplet, \f$\phi_{\rho_3}\f$.
   */
  vector<double> _rho3phase;

  /**
   * maximum weight for the integration of the channel
   */
  vector<double> _maxwgt;

  /**
   * weight for the first  rho multiplet in the integration
   */
  vector<double> _rho1wgt;

  /**
   * weight for the second rho multiplet in the integration
   */
  vector<double> _rho2wgt;

  /**
   * weight for the third  rho multiplet in the integration
   */
  vector<double> _rho3wgt;

  /**
   * mass of the first rho multiplet, \f$M_{\rho_1}\f$.
   */
  vector<double> _rho1mass;

  /**
   * mass of the  second rho multiplet, \f$M_{\rho_2}\f$.
   */
  vector<double> _rho2mass;

  /**
   * mass of the third rho multiplet, \f$M_{\rho_3}\f$.
   */
  vector<double> _rho3mass;

  /**
   * width of the first rho multiplet, \f$\Gamma_{\rho_1}\f$.
   */
  vector<double> _rho1width;

  /**
   * width of the  second rho multiplet, \f$\Gamma_{\rho_2}\f$.
   */
  vector<double> _rho2width;

  /**
   * width of the third rho multiplet, \f$\Gamma_{\rho_3}\f$.
   */
  vector<double> _rho3width;

  /**
   * use the default parameters for the rho masses and widths
   */
  vector<bool> _defaultmass;

  /**
   * constants for the running widths
   */
  vector<vector <double> > _rho0const,_rhocconst;

  /**
   * rho mass parameters
   */
  vector<vector<Energy> > _rhomass;

  /**
   * rho mass parameters
   */
  vector<vector<Energy2> > _rhomass2;

  /**
   * couplings as complex numbers
   */
  vector<vector <Complex> > _ccoupling;

  /**
   * the pion masses
   */
  Energy _mpic,_mpi0;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of VectorMeson3PionDecayer.
 */
template <>
struct BaseClassTrait<Herwig::VectorMeson3PionDecayer,1> {
    /** Typedef of the base class of VectorMeson3PionDecayer. */
  typedef Herwig::VectorMesonDecayerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::VectorMeson3PionDecayer>
  : public ClassTraitsBase<Herwig::VectorMeson3PionDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig++::VectorMeson3PionDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwVMDecay.so"; }

};

}

#include "VectorMeson3PionDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMeson3PionDecayer.tcc"
#endif

#endif /* HERWIG_VectorMeson3PionDecayer_H */
