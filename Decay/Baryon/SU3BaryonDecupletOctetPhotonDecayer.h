// -*- C++ -*-
#ifndef HERWIG_SU3BaryonDecupletOctetPhotonDecayer_H
#define HERWIG_SU3BaryonDecupletOctetPhotonDecayer_H
// This is the declaration of the SU3BaryonDecupletOctetPhotonDecayer class.

#include "Baryon1MesonDecayerBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>SU3BaryonDecupletOctetPhotonDecayer</code> class is designed for the decay
 *  of a deculpet baryon to an octet baryon and a photon using the SU(3) conserving
 *  term in the chiral lagrangian.
 *
 *  The lagrangian is taken to be
 * \f[ \bar{B}_i^b\gamma^\nu\gamma_5\Delta^\mu_{jbc}Q_j\epsilon^{cij}F_{\mu\nu},\f]
 *  where \f$B\f$ is the matrix field for the baryon octet, \f$\Delta\f$ is the field for
 *  the decuplet, \f$Q\f$ is the charge matrix
 * \f[Q  = \left(\begin{array}{ccc}\frac23&0&0\\
 *                                                           0&-\frac13&0\\
 *                                                           0&0&-\frac13
 *     \end{array}\right),\f]
 *  and \f$F^{\mu\nu}\f$ is the electromagnetic field strength. This form is used when 
 *  the multiplets have the same parity and this form without the \f$\gamma_5\f$ if
 *  the multiplets have opposite parity.
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
class SU3BaryonDecupletOctetPhotonDecayer: public Baryon1MesonDecayerBase {

public:

  /**
   * Default constructor.
   */
  SU3BaryonDecupletOctetPhotonDecayer();

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
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac32\f$ and a vector.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A1 The coupling \f$A_1\f$ described above.
   * @param A2 The coupling \f$A_2\f$ described above.
   * @param A3 The coupling \f$A_3\f$ described above.
   * @param B1 The coupling \f$B_1\f$ described above.
   * @param B2 The coupling \f$B_2\f$ described above.
   * @param B3 The coupling \f$B_3\f$ described above.
   */
  virtual void threeHalfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy m2,
					   Complex& A1,Complex& A2,Complex& A3,
					   Complex& B1,Complex& B2,Complex& B3) const;
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
  static ClassDescription<SU3BaryonDecupletOctetPhotonDecayer> initSU3BaryonDecupletOctetPhotonDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  SU3BaryonDecupletOctetPhotonDecayer & operator=(const SU3BaryonDecupletOctetPhotonDecayer &) = delete;

private:

  /**
   * set-up the modes
   */
  void setupModes(unsigned int) const;

private:

  /**
   * the coupling
   */
  InvEnergy C_;

  /**
   * the relative parities of the two baryon multiplets
   */
  bool parity_;

  /**
   * PDG codes for the various octet baryons
   */
  //@{
  /**
   *  The PDG code for the \f$p\f$-like member of the outgoing octet.
   */
  int proton_;

  /**
   *  The PDG code for the \f$n\f$-like member of the outgoing octet.
   */
  int neutron_;

  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the outgoing octet.
   */
  int sigma0_;

  /**
   *  The PDG code for the  \f$\Sigma^+\f$-like member of the outgoing octet.
   */
  int sigmap_;

  /**
   *  The PDG code for the  \f$\Sigma^-\f$-like member of the outgoing octet.
   */
  int sigmam_;

  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the outgoing octet.
   */
  int lambda_;

  /**
   *  The PDG code for the \f$\Xi^0\f$-like member of the outgoing octet.
   */
  int xi0_;

  /**
   *  The PDG code for the \f$\Xi^-\f$-like member of the outgoing octet.
   */
  int xim_;
  //@}

  /**
   * PDG codes for the various decuplet baryons
   */
  //@{
  /**
   *  The PDG code for the \f$\Delta^{++}\f$-like member of the outgoing decuplet.
   */
  int deltapp_;

  /**
   *  The PDG code for the \f$\Delta^{+}\f$-like member of the outgoing decuplet.
   */
  int deltap_;

  /**
   *  The PDG code for the \f$\Delta^{0}\f$-like member of the outgoing decuplet.
   */
  int delta0_;

  /**
   *  The PDG code for the \f$\Delta^{-}\f$-like member of the outgoing decuplet.
   */
  int deltam_;

  /**
   *  The PDG code for the \f$\Sigma^{*+}\f$-like member of the outgoing decuplet.
   */
  int sigmasp_;

  /**
   *  The PDG code for the \f$\Sigma^{*0}\f$-like member of the outgoing decuplet.
   */
  int sigmas0_;

  /**
   *  The PDG code for the \f$\Sigma^{*-}\f$-like member of the outgoing decuplet.
   */
  int sigmasm_;

  /**
   *  The PDG code for the \f$\Omega^-\f$-like member of the outgoing decuplet.
   */
  int omega_;

  /**
   *  The PDG code for the \f$\Xi^{*-}\f$-like member of the outgoing decuplet.
   */
  int xism_;

  /**
   *  The PDG code for the \f$\Xi^{*0}\f$-like member of the outgoing decuplet.
   */
  int xis0_;
  //@}

  /**
   * PDG code for the incoming baryons
   */
  mutable vector<int> incomingB_;

  /**
   * PDG code for the outgoing baryons
   */
  mutable vector<int> outgoingB_;

  /**
   * the maximum weight for the various modes
   */
  vector<double> maxweight_;

  /**
   * The couplings for the different modes.
   */
  mutable vector<InvEnergy> prefactor_;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of SU3BaryonDecupletOctetPhotonDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::SU3BaryonDecupletOctetPhotonDecayer,1> {
    /** Typedef of the base class of SU3BaryonDecupletOctetPhotonDecayer. */
   typedef Herwig::Baryon1MesonDecayerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::SU3BaryonDecupletOctetPhotonDecayer>
  : public ClassTraitsBase<Herwig::SU3BaryonDecupletOctetPhotonDecayer> {
   /** Return the class name. */
  static string className() { return "Herwig::SU3BaryonDecupletOctetPhotonDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwBaryonDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_SU3BaryonDecupletOctetPhotonDecayer_H */
