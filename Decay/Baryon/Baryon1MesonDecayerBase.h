// -*- C++ -*-
#ifndef HERWIG_Baryon1MesonDecayerBase_H
#define HERWIG_Baryon1MesonDecayerBase_H
//
// This is the declaration of the Baryon1MesonDecayerBase class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/PDT/BaryonWidthGenerator.fh"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzSpinor.h"
#include "ThePEG/Helicity/LorentzSpinorBar.h"
#include "ThePEG/Helicity/LorentzRSSpinor.h"
#include "ThePEG/Helicity/LorentzRSSpinorBar.h"
#include "Baryon1MesonDecayerBase.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>Baryon1MesonDecayerBase</code> class is the base class for the decay of
 *  a baryon to another baryon and a pseudoscalar or vector meson.
 *  All the matrix elements involving either a spin-1/2 or spin-3/2 baryons are 
 *  now implemented apart from \f$\frac32\to\frac32+1\f$. The matrix elements
 *  are implemented in a general form with general couplings which must be
 *  supplied in classes inheriting from the one by implementing one of the coupling
 *  members.
 *
 *  - The matrix element for \f$\frac12\to\frac12+0\f$
 *    \f[\mathcal{M} = \bar{u}(p_1)(A+B\gamma_5)u(p_0)\f]
 *  - The matrix element for \f$\frac12\to\frac12+1\f$ 
 *    \f[\mathcal{M} = \bar{u}(p_1)\epsilon^{*\beta}\left[
 *           \gamma_\beta(A_1+B_1\gamma_5)
 *            +p_{0\beta}(A_2+B_2\gamma_5)\right]u(p_0)\f]
 *  - The matrix element for \f$\frac12\to\frac32+0\f$
 *    \f[\bar{u}^\alpha(p_1) p_{0\alpha}\left[A+B\gamma_5\right]u(p_0)\f]
 *  - The matrix element for \f$\frac12\to\frac32+1\f$
 *    \f[\bar{u}^\alpha(p_1)\epsilon^{*\beta}\left[
 *      g_{\alpha\beta}(A_1+B_1\gamma_5)
 *     +p_{0\alpha}(A_2+B_2\gamma_5)
 *     +p_{0\alpha}p_{0\beta}(A_3+B_3\gamma_5)
 *      \right]u(p_0)\f]
 *  - The matrix element for \f$\frac32\to\frac12+0\f$
 *    \f[\bar{u}(p_1) p_{1\alpha}\left[A+B\gamma_5\right]u^\alpha(p_0)\f]
 *  - The matrix element for \f$\frac32\to\frac12+1\f$
 *    \f[\bar{u}(p_1)\epsilon^{*\beta}\left[
 *      g_{\alpha\beta}(A_1+B_1\gamma_5)
 *     +p_{1\alpha}(A_2+B_2\gamma_5)
 *     +p_{1\alpha}p_{0\beta}(A_3+B_3\gamma_5)
 *      \right]u^\alpha(p_0)\f]
 *  - The matrix element for \f$\frac32\to\frac32+0\f$
 *   \f[\bar{u}^\alpha(p_1)\left[(A_1+B_1\gamma_5)g_{\alpha\beta}
 *      +p_{0\alpha}p_{1\beta}(A_2+B_2\gamma_5)\right]u^\beta(p_0)\f]
 *
 * @see DecayIntegrator
 */
class Baryon1MesonDecayerBase: public DecayIntegrator {

  /**
   *  The BaryonWidthGenerator is a friend to get access to the couplings
   */
  friend class BaryonWidthGenerator; 

public:

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * This version uses the generalised couplings to compute the matrix elements
   * given above.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt The option for the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const ParticleVector & decay,MEOption meopt) const;

  /**
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  virtual bool twoBodyMEcode(const DecayMode & dm, int & mecode,
			     double & coupling) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

protected:

  /**
   *  Coupling Members.
   */
  //@{
  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac12\f$ and a scalar.
   * This method must be implemented in any class inheriting from this one
   * which includes \f$\frac12\to\frac12+0\f$ decays. 
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A The coupling \f$A\f$ described above.
   * @param B The coupling \f$B\f$ described above.
   */
  virtual void halfHalfScalarCoupling(int imode, Energy m0, Energy m1, Energy m2,
				      Complex& A,Complex& B) const;

  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac12\f$ and a vector.
   * This method must be implemented in any class inheriting from this one
   * which includes \f$\frac12\to\frac12+1\f$ decays. 
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A1 The coupling \f$A_1\f$ described above.
   * @param A2 The coupling \f$A_2\f$ described above.
   * @param B1 The coupling \f$B_1\f$ described above.
   * @param B2 The coupling \f$B_2\f$ described above.
   */
  virtual void halfHalfVectorCoupling(int imode, Energy m0, Energy m1, Energy m2,
				      Complex& A1,Complex& A2,
				      Complex& B1,Complex& B2) const;

  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac32\f$ and a scalar.
   * This method must be implemented in any class inheriting from this one
   * which includes \f$\frac12\to\frac32+0\f$ or \f$\frac32\to\frac12+0\f$ decays. 
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A The coupling \f$A\f$ described above.
   * @param B The coupling \f$B\f$ described above.
   */
  virtual void halfThreeHalfScalarCoupling(int imode, Energy m0, Energy m1, Energy m2,
					   Complex& A,Complex& B) const;

  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac32\f$ and a vector.
   * This method must be implemented in any class inheriting from this one
   * which includes \f$\frac12\to\frac32+1\f$ or \f$\frac32\to\frac12+1\f$ decays. 
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
  virtual void halfThreeHalfVectorCoupling(int imode, Energy m0, Energy m1, Energy m2,
					   Complex& A1,Complex& A2,Complex& A3,
					   Complex& B1,Complex& B2,Complex& B3) const;


  /**
   * Couplings for spin-\f$\frac32\f$ to spin-\f$\frac12\f$ and a scalar.
   * This method must be implemented in any class inheriting from this one
   * which includes \f$\frac12\to\frac32+0\f$ or \f$\frac32\to\frac12+0\f$ decays. 
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A The coupling \f$A\f$ described above.
   * @param B The coupling \f$B\f$ described above.
   */
  virtual void threeHalfHalfScalarCoupling(int imode, Energy m0, Energy m1, Energy m2,
					   Complex& A,Complex& B) const;

  /**
   * Couplings for spin-\f$\frac32\f$ to spin-\f$\frac12\f$ and a vector.
   * This method must be implemented in any class inheriting from this one
   * which includes \f$\frac12\to\frac32+1\f$ or \f$\frac32\to\frac12+1\f$ decays. 
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
  virtual void threeHalfHalfVectorCoupling(int imode, Energy m0, Energy m1, Energy m2,
					   Complex& A1,Complex& A2,Complex& A3,
					   Complex& B1,Complex& B2,Complex& B3) const;

  /**
   * Couplings for spin-\f$\frac32\f$ to spin-\f$\frac32\f$ and a scalar.
   * This method must be implemented in any class inheriting from this one
   * which includes \f$\frac32\to\frac32+0\f$ decays. 
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

public:

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

private:

  /**
   * Matrix Element Calculation Members
   */
  //@{
  /**
   * Matrix element for spin-\f$\frac12\f$ to spin-\f$\frac12\f$ and a scalar.
   * @param ichan The phase-space channel.
   * @param inpart The decaying particle.
   * @param decay The decay products.
   * @param meopt The option for the matrix element
   * @return The matrix element squared.
   */
  double halfHalfScalar(const int ichan,const Particle & inpart,
			const ParticleVector & decay,MEOption meopt) const;

  /**
   * Matrix element for spin-\f$\frac12\f$ to spin-\f$\frac12\f$ and a vector.
   * @param ichan The phase-space channel.
   * @param inpart The decaying particle.
   * @param decay The decay products.
   * @param meopt The option for the matrix element
   * @return The matrix element squared.
   */
  double halfHalfVector(const int ichan,const Particle & inpart,
			const ParticleVector & decay,MEOption meopt) const;

  /**
   * Matrix element for spin-\f$\frac12\f$ to spin-\f$\frac32\f$ and a scalar.
   * @param ichan The phase-space channel.
   * @param inpart The decaying particle.
   * @param decay The decay products.
   * @param meopt The option for the matrix element
   * @return The matrix element squared.
   */
  double halfThreeHalfScalar(const int ichan,const Particle & inpart,
			     const ParticleVector & decay,MEOption meopt) const;

  /**
   * Matrix element for spin-\f$\frac12\f$ to spin-\f$\frac32\f$ and a vector.
   * @param ichan The phase-space channel.
   * @param inpart The decaying particle.
   * @param decay The decay products.
   * @param meopt The option for the matrix element
   * @return The matrix element squared.
   */
  double halfThreeHalfVector(const int ichan,const Particle & inpart,
			     const ParticleVector  & decay,MEOption meopt) const;

  /**
   * Matrix element for spin-\f$\frac32\f$ to spin-\f$\frac12\f$ and a scalar.
   * @param ichan The phase-space channel.
   * @param inpart The decaying particle.
   * @param decay The decay products.
   * @param meopt The option for the matrix element
   * @return The matrix element squared.
   */
  double threeHalfHalfScalar(const int ichan,const Particle & inpart,
			     const ParticleVector & decay,MEOption meopt) const;

  /**
   * Matrix element for spin-\f$\frac32\f$ to spin-\f$\frac12\f$ and a vector.
   * @param ichan The phase-space channel.
   * @param inpart The decaying particle.
   * @param decay The decay products.
   * @param meopt The option for the matrix element
   * @return The matrix element squared.
   */
  double threeHalfHalfVector(const int ichan,const Particle & inpart,
			     const ParticleVector & decay,MEOption meopt) const;

  /**
   * Matrix element for spin-\f$\frac32\f$ to spin-\f$\frac32\f$ and a scalar.
   * @param ichan The phase-space channel.
   * @param inpart The decaying particle.
   * @param decay The decay products.
   * @param meopt The option for the matrix element
   * @return The matrix element squared.
   */
  double threeHalfThreeHalfScalar(const int ichan,const Particle & inpart,
				  const ParticleVector & decay,MEOption meopt) const;
  //@}

private:

  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractNoPIOClassDescription<Baryon1MesonDecayerBase> 
  initBaryon1MesonDecayerBase;

  /**
   * Private and non-existent assignment operator.
   */
  Baryon1MesonDecayerBase & operator=(const Baryon1MesonDecayerBase &) = delete;

private:

  /**
   *  Spin density matrx
   */
  mutable RhoDMatrix _rho;

  /**
   *  Spin-\f$\frac12\f$ spinor
   */
  mutable vector<Helicity::LorentzSpinor<SqrtEnergy> >      _inHalf;

  /**
   *  Spin-\f$\frac12\f$ barred spinor
   */
  mutable vector<Helicity::LorentzSpinorBar<SqrtEnergy> >   _inHalfBar;

  /**
   *  Spin-\f$\frac32\f$ spinor
   */
  mutable vector<Helicity::LorentzRSSpinor<SqrtEnergy> >    _inThreeHalf;

  /**
   *  Spin-\f$\frac32\f$ barred spinor
   */
  mutable vector<Helicity::LorentzRSSpinorBar<SqrtEnergy> > _inThreeHalfBar;

  /**
   *  Polarization vector
   */
  mutable vector<Helicity::LorentzPolarizationVector> _inVec;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of Baryon1MesonDecayerBase.
 */
template <>
 struct BaseClassTrait<Herwig::Baryon1MesonDecayerBase,1> {
    /** Typedef of the base class of Baryon1MesonDecayerBase. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Baryon1MesonDecayerBase>
  : public ClassTraitsBase<Herwig::Baryon1MesonDecayerBase> {
  /** Return the class name. */
  static string className() { return "Herwig::Baryon1MesonDecayerBase";}
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwBaryonDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_Baryon1MesonDecayerBase_H */
