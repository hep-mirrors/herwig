// -*- C++ -*-
#ifndef HERWIG_SU3BaryonOctetDecupletScalarDecayer_H
#define HERWIG_SU3BaryonOctetDecupletScalarDecayer_H
// This is the declaration of the SU3BaryonOctetDecupletScalarDecayer class.

#include "Baryon1MesonDecayerBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 * 
 *  The <code>SU3BaryonOctetDecupletScalarDecayer</code> class is designed for the decay
 *  of an SU(3) octet baryon to a deculpet baryon and a pseudoscalar meson from the
 *  lightest multiplet.
 *
 *  The coupling is taken to have the form
 *  \f[\frac{C}{f_\pi}\left[
 *  \bar{\Delta}^{\mu,abc}p_{\phi,\mu}\phi^i_aB^j_b\epsilon_{cij}
 * -\bar{B}^b_ip_{\phi,\mu}\phi^a_j\Delta^\mu_{abc}\epsilon^{cij}
 *  \right],\f]
 *  where \f$\Delta^\mu_{abc}\f$ is the decuplet field, \f$B^j_b\f$ is the
 *  octet field,
 *  \f$\phi^i_a\f$ is the pseudoscalar field, \f$f_\pi\f$ is the pion decay constant
 *  and \f$C\f$ is the coupling for the decay.
 *  This form is used for the \f$\frac12^+\f$ octet mesons and the same form with
 *  \f$\gamma_5\f$ for the \f$\frac12^-\f$ octet.
 *
 *  The above form is used for the decay \f$\frac12\to\frac32+0\f$ whereas for the
 *  decay of a spin-\f$\frac32\f$ octet we use
 *  \f[\frac{C}{f_\pi}\left[
 *  \bar{\Delta}^{\mu,abc}p_\phi\!\!\!\!\!\!\!\!\not\,\,\,\,\,\phi^i_aB^j_{b\mu}\epsilon_{cij}
 * -\bar{B}^b_{i\mu}p_\phi\!\!\!\!\!\!\!\!\not\,\,\,\,\,\phi^a_j\Delta^\mu_{abc}\epsilon^{cij}
 *  \right]\f]
 *  where the decuplet and octet have the same parity and the same form with the 
 *  addition of \f$\gamma_5\f$ if they have opposite parity.
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
class SU3BaryonOctetDecupletScalarDecayer: public Baryon1MesonDecayerBase {

public:

  /**
   * Default constructor.
   */
  SU3BaryonOctetDecupletScalarDecayer();

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
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac32\f$ and a scalar.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A The coupling \f$A\f$ described above.
   * @param B The coupling \f$B\f$ described above.
   */
  virtual void halfThreeHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
					   Complex& A,Complex& B) const;

  /**
   * Couplings for spin-\f$\frac32\f$ to spin-\f$\frac32\f$ and a scalar.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A1 The coupling \f$A_1\f$ described above.
   * @param A2 The coupling \f$A_2\f$ described above.
   * @param B1 The coupling \f$B_1\f$ described above.
   * @param B2 The coupling \f$B_2\f$ described above.
   */
  virtual void threeHalfThreeHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
						Complex& A1,Complex& A2,
						Complex& B1,Complex& B2) const;
  //@}

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
  static ClassDescription<SU3BaryonOctetDecupletScalarDecayer> initSU3BaryonOctetDecupletScalarDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  SU3BaryonOctetDecupletScalarDecayer & operator=(const SU3BaryonOctetDecupletScalarDecayer &);

private:

  /**
   * set-up the modes
   */
  void setupModes(unsigned int) const;

private:

  /**
   * the couplings
   */
  double _c;

  /**
   * the relative parities of the two baryon multiplets
   */
  bool _parity;

  /**
   * the pion decay constant
   */
  Energy _fpi;

  /**
   * PDG codes for the various octet baryons
   */
  //@{
  /**
   *  The PDG code for the \f$p\f$-like member of the incoming octet.
   */
  int _proton;

  /**
   *  The PDG code for the \f$n\f$-like member of the incoming octet.
   */
  int _neutron;

  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the incoming octet.
   */
  int _sigma0;

  /**
   *  The PDG code for the  \f$\Sigma^+\f$-like member of the incoming octet.
   */
  int _sigmap;

  /**
   *  The PDG code for the  \f$\Sigma^-\f$-like member of the incoming octet.
   */
  int _sigmam;

  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the incoming octet.
   */
  int _lambda;

  /**
   *  The PDG code for the \f$\Xi^0\f$-like member of the incoming octet.
   */
  int _xi0;

  /**
   *  The PDG code for the \f$\Xi^-\f$-like member of the incoming octet.
   */
  int _xim;
  //@}

  /**
   * PDG codes for the various decuplet baryons
   */
  //@{
  /**
   *  The PDG code for the \f$\Delta^{++}\f$-like member of the outgoing decuplet.
   */
  int _deltapp;

  /**
   *  The PDG code for the \f$\Delta^{+}\f$-like member of the outgoing decuplet.
   */
  int _deltap;

  /**
   *  The PDG code for the \f$\Delta^{0}\f$-like member of the outgoing decuplet.
   */
  int _delta0;

  /**
   *  The PDG code for the \f$\Delta^{-}\f$-like member of the outgoing decuplet.
   */
  int _deltam;

  /**
   *  The PDG code for the \f$\Sigma^{*+}\f$-like member of the outgoing decuplet.
   */
  int _sigmasp;

  /**
   *  The PDG code for the \f$\Sigma^{*0}\f$-like member of the outgoing decuplet.
   */
  int _sigmas0;

  /**
   *  The PDG code for the \f$\Sigma^{*-}\f$-like member of the outgoing decuplet.
   */
  int _sigmasm;

  /**
   *  The PDG code for the \f$\Omega^-\f$-like member of the outgoing decuplet.
   */
  int _omega;

  /**
   *  The PDG code for the \f$\Xi^{*-}\f$-like member of the outgoing decuplet.
   */
  int _xism;

  /**
   *  The PDG code for the \f$\Xi^{*0}\f$-like member of the outgoing decuplet.
   */
  int _xis0;
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
   * The couplings for the different modes.
   */
  mutable vector<InvEnergy> _prefactor;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of SU3BaryonOctetDecupletScalarDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::SU3BaryonOctetDecupletScalarDecayer,1> {
    /** Typedef of the base class of SU3BaryonOctetDecupletScalarDecayer. */
   typedef Herwig::Baryon1MesonDecayerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::SU3BaryonOctetDecupletScalarDecayer>
    /** Typedef of the base class of SU3BaryonOctetDecupletScalarDecayer. */
  : public ClassTraitsBase<Herwig::SU3BaryonOctetDecupletScalarDecayer> {
   /** Return the class name.*/
  static string className() { return "Herwig::SU3BaryonOctetDecupletScalarDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwBaryonDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_SU3BaryonOctetDecupletScalarDecayer_H */
