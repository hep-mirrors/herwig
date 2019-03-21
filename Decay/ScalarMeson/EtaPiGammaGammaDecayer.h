// -*- C++ -*-
//
// EtaPiGammaGammaDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_EtaPiGammaGammaDecayer_H
#define HERWIG_EtaPiGammaGammaDecayer_H
// This is the declaration of the EtaPiGammaGammaDecayer class.

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>EtaPiGammaGammaDecayer</code> class implements a VMD model
 * matrix element for \f$\eta,\eta'\to \pi^0 \gamma \gamma\f$ taken from
 * hep-ph/0112150.
 *
 *  The matrix element is given by
 *  \f[ \mathcal{M} = 
 *    D(s,t,u)\left[ \epsilon_1\cdot\epsilon_2 q_1\cdot q_2
 *                  -\epsilon_1\cdot q_2\epsilon_2\cdot q_1\right]
 *   -E(s,t,u)\left[-\epsilon_1\cdot\epsilon_2p\cdot q_1p\cdot q_2
 *                  -\epsilon_1\cdot p\epsilon_2\cdot p q_1\cdot q_2
 *                  +\epsilon_1\cdot q_2 \epsilon_2\cdot  p  p\cdot q_1
 *                  +\epsilon_1\cdot  p  \epsilon_2\cdot q_1 p\cdot q_2
 *\right],
 *  \f]
 *  where \f$q_{1,2}\f$ are the momenta of the photons, \f$\epsilon_{1,2}\f$ are
 *  the polarization vectors of the photons, \f$p\f$ is the momentum of the decaying
 *  meson and
 * \f[D(s,t,u) = \frac{2\sqrt{3}}{9}g^2_{\omega\rho\pi}
 *              \left(\frac{2eF^2_{\pi}g}{m^2_V}\right)^2\times\left(\frac{F_\pi}{F_8}\cos\theta\mp\sqrt{2}\frac{F_\pi}{F_0}\sin\theta\right)\times\left[\frac{p\cdot q_2-m^2_\eta}{m^2_V-t}+\frac{p\cdot q_1-m^2_\eta}{m^2_V-u}\right]\f]
 * \f[E(s,t,u) =-\frac{2\sqrt{3}}{9}g^2_{\omega\rho\pi}\left(\frac{2eF^2_{\pi}g}{m^2_V}\right)^2\times\left(\frac{F_\pi}{F_8}\cos\theta\mp\sqrt{2}\frac{F_\pi}{F_0}\sin\theta\right)\times\left[\frac1{m^2_V-t}+\frac1{m^2_V-u}\right].\f]
 *  the \f$-\f$ sign corresponds to the \f$\eta\f$ decay and the \f$+\f$ to the \f$\eta'\f$ decay.
 *  Here
 * - \f$g_{\omega\rho\pi}\f$ is the coupling of the
 *               \f$\omega\f$ to the \f$\rho\f$ and a pion.
 * - \f$e\f$ is the electric charge.
 * - \f$F_{\pi}\f$ is the pion decay constant.
 * - \f$g\f$ is the conversion factor for a \f$\rho\f$ into a photon.
 * - \f$m_V\f$ is the mass of the vector meson, in this case we use the \f$\rho\f$.
 * - \f$F_0\f$ is the singlet decay constant.
 * - \f$F_8\f$ is the octet decay constant.
 * - \f$\theta\f$ is the octet-singlet mixing angle
 * - \f$m_\eta\f$ is the mass of the decay meson.
 *
 *  In practice we use a slightly modified form by including a running width term to 
 * include the \f$\eta'\f$ decay as well as the \f$\eta\f$ decay.
 *
 * @see DecayIntegrator
 * 
 */
class EtaPiGammaGammaDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  EtaPiGammaGammaDecayer();
  
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
  virtual double threeBodyMatrixElement(const int imode, const Energy2 q2,
					const  Energy2 s3, const Energy2 s2,
					const Energy2 s1, const Energy m1,
					const Energy m2,const Energy m3) const;

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
  static ClassDescription<EtaPiGammaGammaDecayer> initEtaPiGammaGammaDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  EtaPiGammaGammaDecayer & operator=(const EtaPiGammaGammaDecayer &) = delete;

private:

  /**
   * The coupling \f$g_{\omega\rho\pi}\f$ of the \f$\rho\f$ to \f$\omega\f$ and
   * a \f$\pi\f$.
   */
  InvEnergy _grhoomega;

  /**
   * The pion decay constant, \f$F_{\pi}\f$
   */
  Energy _fpi;

  /**
   * The mass of the \f$\rho\f$.
   */
  Energy _rhomass;

  /**
   * The width of the \f$\rho\f$. 
   */
  Energy _rhowidth;

  /**
   * The coupling for the conversion of a rho to a photon, \f$g\f$
   */
  double _grho;

  /**
   * The mass of the pion
   */
  Energy _mpi;

  /**
   * Constant for the running \f$\rho\f$ width calculation.
   */
  double _rhoconst;

  /**
   * Use local values of the \f$\rho\f$ mass and width.
   */
  bool _localparameters;

  /**
   * Ratios of the decay constants \f$F_8/F_\pi\f$.
   */
  double _ratiofpif8;

  /**
   * Ratios of the decay constants \f$F_0/F_\pi\f$.
   */
  double _ratiofpif0;

  /**
   * the mixing angle, \f$\theta\f$.
   */
  double _theta;

  /**
   * the maximum weights for the \f$\eta\f$ decay.
   */
  double _etamax;

  /**
   * the maximum weights for the \f$\eta'\f$ decay.
   */
  double _etapmax;

  /**
   * The prefactor for the \f$D(s,t,u)\f$ function.
   */
  vector<InvEnergy2> _dconst;

  /**
   * The prefactor for the \f$E(s,t,u)\f$ function.
   */
  vector<InvEnergy2> _econst;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization vectors for the photons
   */
  mutable vector<Helicity::LorentzPolarizationVector> _vectors[2];
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of EtaPiGammaGammaDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::EtaPiGammaGammaDecayer,1> {
    /** Typedef of the base class of EtaPiGammaGammaDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::EtaPiGammaGammaDecayer>
  : public ClassTraitsBase<Herwig::EtaPiGammaGammaDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig::EtaPiGammaGammaDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwSMDecay.so"; }
  
};

/** @endcond */
  
}

#endif /* HERWIG_EtaPiGammaGammaDecayer_H */
