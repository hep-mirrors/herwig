// -*- C++ -*-
#ifndef THEPEG_SU3BaryonSingletOctetPhotonDecayer_H
#define THEPEG_SU3BaryonSingletOctetPhotonDecayer_H
// This is the declaration of the SU3BaryonSingletOctetPhotonDecayer class.

#include "Baryon1MesonDecayerBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>SU3BaryonSingletOctetPhotonDecayer</code> class performs the
 *  radiative decay  of an \f$SU(3)\f$ singlet baryon to an \f$SU(3)\f$ octet baryon.
 *
 * The Lagrangian is taken to be
 * \f[iC\left[ {\rm tr}\left(\bar{S}\sigma_{\mu\nu}f^{\mu\nu}_+,B\right)
 *            +{\rm tr}\left(\bar{B}\sigma_{\mu\nu}f^{\mu\nu}_+,S\right)
 *      \right]
 * \f]
 *  where \f$B\f$ is the matrix field for the ground state baryon multiplet, \f$S\f$
 *  is the singlet field and \f$f^{\mu\nu}_+\f$ is the chiral
 * field strength tensor for the electromagentic field given by
 * \f[ f^{\mu\nu}_+ = Q F^{\mu\nu} = \left(\begin{array}{ccc}\frac23&0&0\\
 *                                                           0&-\frac13&0\\
 *                                                           0&0&-\frac13
 *     \end{array}\right) F^{\mu\nu},\f]
 * where \f$F^{\mu\nu}\f$ is the electromagentic field strength tensor.
 * This form is used for the case where both the singlet and the baryon multiplet
 * have the same parity and
 * an additional \f$\gamma_5\f$ is added for the case where the multiplets have 
 * opposite parity.
 *
 * For the decay of spin-\f$\frac32\f$ baryons we use the form
 * \f[
 *  ir_d\left[ {\rm tr}\left(\bar{S}^\mu\gamma_\nu f^{\mu\nu}_+,B\right)
 *            +{\rm tr}\left(\bar{B}\gamma_\nu f^{\mu\nu}_+,S^\mu\right)
 *            \right],\f]
 *  where \f$S^\mu\f$ is the matrix field for the excited baryon multiplet.
 *  This form is used when the baryon's have the same parity and this form with
 *  an additional \f$\gamma_5\f$ when they have opposite parity.
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
 */
class SU3BaryonSingletOctetPhotonDecayer: public Baryon1MesonDecayerBase {

public:

  /**
   * Default constructor.
   */
  SU3BaryonSingletOctetPhotonDecayer();

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
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac12\f$ and a vector.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A1 The coupling \f$A_1\f$ described above.
   * @param A2 The coupling \f$A_2\f$ described above.
   * @param B1 The coupling \f$B_1\f$ described above.
   * @param B2 The coupling \f$B_2\f$ described above.
   */
  virtual void halfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy m2,
				      Complex& A1,Complex& A2,
				      Complex& B1,Complex& B2) const;

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
  static ClassDescription<SU3BaryonSingletOctetPhotonDecayer> initSU3BaryonSingletOctetPhotonDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  SU3BaryonSingletOctetPhotonDecayer & operator=(const SU3BaryonSingletOctetPhotonDecayer &) = delete;

private:

  /**
   * set-up the modes
   */
  void setupModes(unsigned int) const;

private:

  /**
   * the coupling
   */
  InvEnergy _c;

  /**
   * the relative parities of the two baryon multiplets
   */
  bool _parity;

  /**
   * PDG codes for the lower lying baryon octet baryons
   */
  //@{
  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the outgoing octet.
   */
  int _sigma0;

  /**
   *  The PDG code for the \f$\Sigma^0\f$-like member of the outgoing octet.
   */
  int _lambda;

  /**
   * PDG code for the various excited baryon
   */
  int _elambda;

  /**
   * outgoing baryon for the various modes
   */
  mutable vector<int> _outgoingB;

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
 * base class of SU3BaryonSingletOctetPhotonDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::SU3BaryonSingletOctetPhotonDecayer,1> {
    /** Typedef of the base class of  SU3BaryonSingletOctetPhotonDecayer. */
   typedef Herwig::Baryon1MesonDecayerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::SU3BaryonSingletOctetPhotonDecayer>
  : public ClassTraitsBase<Herwig::SU3BaryonSingletOctetPhotonDecayer> {
   /** Return the class name.*/
  static string className() { return "Herwig::SU3BaryonSingletOctetPhotonDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwBaryonDecay.so"; }

};

/** @endcond */

}

#endif /* THEPEG_SU3BaryonSingletOctetPhotonDecayer_H */
