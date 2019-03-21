// -*- C++ -*-
#ifndef HERWIG_SU3BaryonOctetOctetScalarDecayer_H
#define HERWIG_SU3BaryonOctetOctetScalarDecayer_H
// This is the declaration of the SU3BaryonOctetOctetScalarDecayer class.

#include "Baryon1MesonDecayerBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>SU3BaryonOctetOctetScalarDecayer</code> class is a simple decayer for 
 * the strong decay of the excited baryon \f$SU(3)\f$ octets to lower lying
 * baryon octets and a pseudoscalar meson from the lightest 
 * multiplet, i.e. \f$\pi\f$, \f$K\f$, \f$\eta\f$, etc..
 *
 * The interaction Lagrangian is taken to have the form
 *
 * \f[-\frac{s_d}{4\sqrt{2}f_\pi} \left[{\rm tr} (\bar{R}p_\phi\!\!\!\!\!\!\!\!\not\,\,\,\,\,
 *     \{\phi,B\})\right]
 * -\frac{s_f}{4\sqrt{2}f_\pi} \left[{\rm tr}    (\bar{R}p_\phi\!\!\!\!\!\!\!\!\not\,\,\,\,\,
 *      [\phi,B])\right],\f]
 *  where \f$R\f$ is the matrix field for the excited resonance, \f$\phi\f$ is the matrix
 *  field for the pseudoscalar mesons, \f$f_\pi\f$ is the pion decay constant
 *  and \f$B\f$ is the matrix field for the ground
 *  state baryon octet for the decay of a spin\f$\frac12\f$ multiplet.
 *
 *  The above form is used for \f$\frac12^-\f$ excited baryon resonances
 *  and the above form with
 *  an additional \f$\gamma_5\f$ is used for excited \f$\frac12^+\f$ baryon resonances.

 *  For the decay of a spin-\f$\frac32\f$ multiplet we use the form 
 *  
 * \f[-\frac{s_d}{4\sqrt{2}f_\pi} \left[{\rm tr} (\bar{R}^\mu p_{\phi,\mu}\{\phi,B\})\right]
 *    -\frac{s_f}{4\sqrt{2}f_\pi} \left[{\rm tr} (\bar{R}^\mu p_{\phi,\mu} [\phi,B])\right],\f]
 *  where \f$R^\mu\f$ is the matrix field for the excited resonance. This form is used
 *  for \f$\frac32^+\f$ excited baryon resonances and this form with an additional
 *  \f$\gamma_5\f$ for excited \f$\frac32^-\f$ resonances.
 *
 *  The decay to spin-\f$\frac32\f$ resonances in not yet implemented but can
 *  be included if necessary.
 *
 *  This is one of a number of decayers based on \f$SU(3)\f$ symmetry which are
 *  intended for the decay of excited baryons.
 *
 * @see Baryon1MesonDecayerBase
 * @see SU3BaryonOctetOctetPhotonDecayer
 * @see SU3BaryonDecupletOctetPhotonDecayer
 * @see SU3BaryonDecupletOctetScalarDecayer
 * @see SU3BaryonSingletOctetPhotonDecayer
 * @see SU3BaryonSingletOctetScalarDecayer
 * @see SU3BaryonOctetDecupletScalarDecayer
 * @see SU3BaryonOctetOctetScalarDecayer
 * 
 */
class SU3BaryonOctetOctetScalarDecayer: public Baryon1MesonDecayerBase {

public:

  /**
   * Default constructor.
   */
  SU3BaryonOctetOctetScalarDecayer();

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

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

  /**
   *  Coupling Members.
   */
  //@{
  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac12\f$ and a scalar.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A The coupling \f$A\f$ described above.
   * @param B The coupling \f$B\f$ described above.
   */
  virtual void halfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
				      Complex& A,Complex& B) const;

  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac32\f$ and a scalar.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A The coupling \f$A\f$ described above.
   * @param B The coupling \f$B\f$ described above.
   */
  virtual void threeHalfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
					   Complex& A,Complex& B) const;
  //@}

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
   * Private and non-existent assignment operator.
   */
  SU3BaryonOctetOctetScalarDecayer & operator=(const SU3BaryonOctetOctetScalarDecayer &) = delete;

private:

  /**
   * set-up the modes
   */
  void setupModes(unsigned int) const;

private:

  /**
   * The couplings of comutator term.
   */
  double _sf;

  /**
   * The couplings of anticomutator term.
   */
  double _sd;

  /**
   * the relative parities of the two baryon multiplets
   */
  bool _parity;

  /**
   * the pion decay constant
   */
  Energy _fpi;

  /**
   * PDG codes for the lower lying baryon octet baryons
   */
  //@{
  /**
   *  The PDG code for the \f$p\f$-like member of the outgoing octet.
   */
  int _proton;

  /**
   *  The PDG code for the \f$n\f$-like member of the outgoing octet.
   */
  int _neutron;

  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the outgoing octet.
   */
  int _sigma0;

  /**
   *  The PDG code for the  \f$\Sigma^+\f$-like member of the outgoing octet.
   */
  int _sigmap;

  /**
   *  The PDG code for the  \f$\Sigma^-\f$-like member of the outgoing octet.
   */
  int _sigmam;

  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the outgoing octet.
   */
  int _lambda;

  /**
   *  The PDG code for the \f$\Xi^0\f$-like member of the outgoing octet.
   */
  int _xi0;

  /**
   *  The PDG code for the \f$\Xi^-\f$-like member of the outgoing octet.
   */
  int _xim;
  //@}


  /**
   * PDG codes for the higher  baryon octet baryons
   */
  //@{
  /**
   *  The PDG code for the \f$p\f$-like member of the outgoing octet.
   */
  int _eproton;

  /**
   *  The PDG code for the \f$n\f$-like member of the outgoing octet.
   */
  int _eneutron;

  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the outgoing octet.
   */
  int _esigma0;

  /**
   *  The PDG code for the  \f$\Sigma^+\f$-like member of the outgoing octet.
   */
  int _esigmap;

  /**
   *  The PDG code for the  \f$\Sigma^-\f$-like member of the outgoing octet.
   */
  int _esigmam;

  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the outgoing octet.
   */
  int _elambda;

  /**
   *  The PDG code for the \f$\Xi^0\f$-like member of the outgoing octet.
   */
  int _exi0;

  /**
   *  The PDG code for the \f$\Xi^-\f$-like member of the outgoing octet.
   */
  int _exim;
  //@}

  /**
   * PDG code for the incoming baryons
   */
  mutable vector<int> _incomingB;

  /**
   * PDG code for the outgoing baryons
   */
  mutable vector<int> _outgoingB;

  /**
   * PDG code for the outgoing mesons
   */
  mutable vector<int> _outgoingM;

  /**
   * the maximum weight for the various modes
   */
  vector<double> _maxweight;

  /**
   * The couplings for the different modes
   */
  mutable vector<InvEnergy> _prefactor;

};

}


#endif /* HERWIG_SU3BaryonOctetOctetScalarDecayer_H */
