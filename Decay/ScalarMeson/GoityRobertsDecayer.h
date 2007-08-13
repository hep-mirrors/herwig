// -*- C++ -*-
#ifndef HERWIG_GoityRobertsDecayer_H
#define HERWIG_GoityRobertsDecayer_H
//
// This is the declaration of the GoityRobertsDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/WeakCurrents/LeptonNeutrinoCurrent.h"
#include "GoityRobertsDecayer.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The GoityRobertsDecayer class implements the model of PRD51, 3459 (1995) for
 * the semileptonic decay of B mesons to \f$D^{(*)}\pi\ell\nu\f$.
 *
 * @see DecayIntegrator
 * @see LeptonNeutrinoCurrent 
 * @see \ref GoityRobertsDecayerInterfaces "The interfaces"
 * defined for GoityRobertsDecayer.
 */
class GoityRobertsDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  GoityRobertsDecayer();
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const PDVector & children) const;
  
  /**
   * Check if this decayer can perfom the decay for a particular mode.
   * Uses the modeNumber member but can be overridden
   * @param parent The decaying particle
   * @param children The decay products
   */
  inline virtual bool accept(tcPDPtr parent, const PDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * This function combines the current and the form factor to give the matrix
   * element.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(bool vertex, const int ichan, const Particle & part,
		     const ParticleVector & decay) const;

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<GoityRobertsDecayer> initGoityRobertsDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GoityRobertsDecayer & operator=(const GoityRobertsDecayer &);

private:

  /**
   * Calculate the Isgur-Wise functions for the decays.
   * @param omega The velocity transfer, \f$\omega\f$.
   * @param xi    The \f$\xi(\omega)\f$ form factor.
   * @param xi1   The \f$\xi^(1)(\omega)\f$ form factor.
   * @param rho1  The \f$\rho_1(\omega)\f$ form factor.
   * @param rho2  The \f$\rho_2(\omega)\f$ form factor.
   */
  inline void calculateFormFactors(double omega,double & xi,double & xi1,
				   double & rho1,double & rho2) const;

private:

  /**
   * The leptonic current
   */
  LeptonNeutrinoCurrentPtr _current;

  /**
   * the fermi constant
   */
  InvEnergy2 _GF;

  /**
   *  Include the \f$D^*\f$ in the \f$B\to D\pi\f$ decay.
   */
  bool _includeDstar;

  /**
   *  Include the non-\f$D^*\f$ excited states in the \f$B\to D^{(*)}\pi\f$ decay.
   */
  bool _includeDstarstar;

  /**
   * Wavefunction \f$\beta\f$ parameters for the 1S mesons.
   */
  Energy _beta1S;

  /**
   * Wavefunction \f$\beta\f$ parameter for the 2S mesons.
   */
  Energy _beta2S;

  /**
   * Wavefunction \f$\beta\f$ parameter for the 1P \f$(0^+,1^+)\f$ states.
   */
  Energy _beta1P;

  /**
   * Wavefunction \f$\beta\f$ parameter for the 1D \f$(1^-,2^-)\f$ states.
   */
  Energy _beta1D;

  /**
   *  The pion decay constant, \f$f_\pi\f$.
   */
  Energy _fpi;

  /**
   *  The mass difference for the 2S mesons.
   */
  Energy _deltaM2S;

  /**
   * The mass difference for the 1P mesons.
   */
  Energy _deltaM1P;

  /**
   * The mass difference for the 1D mesons
   */
  Energy _deltaM1D;

  /**
   * The width for the 2S mesons
   */ 
  Energy _gamma2S;

  /**
   * The width for the 1P mesons
   */ 
  Energy _gamma1P;

  /**
   * The width for the 1D mesons
   */ 
  Energy _gamma1D;

  /**
   * The \f$\bar{\Lambda}\f$ parameter for the form factors.
   */
  Energy _lambdabar;

  /**
   * The \f$g\f$ coupling for the decays.
   */
  double _g;

  /**
   * The \f$\alpha_1\f$ coupling for the decays.
   */
  double _alpha1;

  /**
   * The \f$\alpha_2\f$ coupling for the decays
   */
  double _alpha2;

  /**
   * The \f$\alpha_3\f$ coupling for the decays
   */
  double _alpha3;

  /**
   * location of the weights
   */
  vector<int> _wgtloc;

  /**
   * the maximum weights
   */
  vector<double> _wgtmax;

  /**
   *  The weights for the different channels
   */
  vector<double> _weights;

  /**
   *  Parameters for the masses and mass differences
   */
  //@{
  /**
   *  Mass difference between ground and excited B mesons
   */
  Energy _deltaMb;

  /**
   *  Width of the excited B mesons
   */
  Energy _gammaB;

  /**
   *  Width of the \f$D^{*0}\f$
   */
  Energy _gammaD0;

  /**
   *  Width of the \f$D^{*+}\f$
   */
  Energy _gammaDp;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GoityRobertsDecayer. */
template <>
struct BaseClassTrait<Herwig::GoityRobertsDecayer,1> {
  /** Typedef of the first base class of GoityRobertsDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GoityRobertsDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GoityRobertsDecayer>
  : public ClassTraitsBase<Herwig::GoityRobertsDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::GoityRobertsDecayer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * GoityRobertsDecayer is implemented. It may also include several, space-separated,
   * libraries if the class GoityRobertsDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSMDecay.so"; }
};

/** @endcond */

}

#include "GoityRobertsDecayer.icc"

#endif /* HERWIG_GoityRobertsDecayer_H */
