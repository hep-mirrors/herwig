// -*- C++ -*-
#ifndef HERWIG_FFDipole_H
#define HERWIG_FFDipole_H
//
// This is the declaration of the FFDipole class.
//

#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "Herwig++/Utilities/Math.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "ThePEG/Interface/Interfaced.h"
#include "FFDipole.fh"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Constants::pi;

/** \ingroup Decay
 *
 * The FFDipole class generates radiation from a final-final dipole for
 * the generation of photons in decay by the YODA algorithm.
 * 
 * @see YODA
 */
class FFDipole: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline FFDipole();
  //@}

public:

  /**
   * Member to generate the photons from the dipole
   * @param p The decaying particle
   * @param children The decay products
   * @return The decay products with additional radiation
   */
  virtual ParticleVector generatePhotons(const Particle & p,
					 ParticleVector children);

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

  /**
   * Return the photon multiplicity according to a Poissonian Distribution
   * with the supplied average
   * @param average The average
   * @return A value from tghe poisson distribution
   */
  inline int poisson(double average);

  /**
   * Generate the momentum of a photon 
   * @param beta1 The velocity, \f$\beta_1\f$, of the first charged particle
   * @param ombeta1 One minus the velocity, \f$1-\beta_1\f$, of the first 
   * charged particle which is supplied for numerical stability
   * @param beta2 The velocity, \f$\beta_2\f$, of the second charged particle
   * @param ombeta2 One minus the velocity, \f$1-\beta_2\f$, of the 
   * second charged particle which is supplied for numerical stability
   * @return The contribution to the dipole weight
   */
  double photon(double beta1,double ombeta1, double beta2, double ombeta2);

  /**
   * Calculate the exact weight for the dipole.
   * @param beta1 Velocity of the first charged particle, \f$\beta_1\f$
   * @param beta2 Velocity of the second charged particle, \f$\beta_2\f$.
   * @param ombeta1 One minus the velocity of the first particle,  \f$1-\beta_1\f$
   * @param ombeta2 One minus the velocity of the second particle,  \f$1-\beta_2\f$
   * @param iphot The number of the photon for which the weight is required
   * @return The weight
   */
  inline double exactDipoleWeight(double beta1,double ombeta1,
				  double beta2,double ombeta2,unsigned int iphot);

  /**
   * Jacobian factor for the weight
   */
  inline double jacobianWeight();

  /**
   * Matrix element weight
   */
  double meWeight(ParticleVector children);

  /**
   * Member which generates the photons
   * @param boost Boost vector to take the particles produced back from
   * the decaying particle's rest frame to the lab
   * @param children The decay products
   */
  double makePhotons(Boost boost,ParticleVector children);

  /**
   *  Boost all the momenta from the dipole rest frame via the parent rest frame
   * to the lab
   * @param boost The boost vector from the rest frame to the lab
   * @return Whether or not it suceeded
   */
  bool boostMomenta(Boost boost);

  /**
   *  Remove any photons which fail the energy cuts
   * @return Number of photons removed
   */
  unsigned int removePhotons();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<FFDipole> initFFDipole;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FFDipole & operator=(const FFDipole &);

private:

  /**
   *  The minimum photon energy in the boosted frame
   */
  Energy _emin;

  /**
   *  The minimum photon energy in the rest frame
   */
  Energy _eminrest;

  /**
   *  The minimum photon energy in the lab  frame
   */
  Energy _eminlab;

  /**
   *  The maximum photon energy
   */
  Energy _emax;

  /**
   *  Photon multiplicity being generated
   */
  unsigned int _multiplicity;

  /**
   * Maximum number of photons to generate
   */
  unsigned int _nphotonmax;

  /**
   *  Masses of the particles involved
   */
  vector<Energy> _m;

  /**
   *  Produce of the particles charges
   */
  double _charge;

  /**
   *   Momenta of the particles in the dipole rest frame
   */
  //@{
  /**
   *  Momenta of the charged particles in the dipole rest frame before radiation
   */
  vector<Lorentz5Momentum> _qdrf;

  /**   *  Momenta of the charged particles in the dipole rest frame after radiation
   */
  vector<Lorentz5Momentum> _qnewdrf;

  /**
   *  Momenta of the photons in the dipole rest frame
   */
  vector<Lorentz5Momentum> _ldrf;

  /**
   * Total momentum of the photons in the dipole rest frame
   */
  Lorentz5Momentum _bigLdrf;
  //@}

  /**
   *  Momentum of the particles in the parent's rest frame
   */
  //@{
  /**
   *  Momenta of the charged particles in the parent's rest frame before radiation
   */
  vector<Lorentz5Momentum> _qprf;

  /**
   *  Momenta of the charged particles in the parent's rest frame after radiation
   */
  vector<Lorentz5Momentum> _qnewprf;

  /**
   *  Momenta of the photons in the parent rest frame
   */
  vector<Lorentz5Momentum> _lprf;

  /**
   * Total momentum of the photons in the parent rest frame
   */
  Lorentz5Momentum _bigLprf;
  //@}

  /**
   *  Momentum of the particles in the lab frame
   */
  //@{
  /**
   *  Momenta of the charged particles in the lab frame before radiation
   */
  vector<Lorentz5Momentum> _qlab;
  
  /**
   *  Momenta of the charged particles in the lab frame after  radiation
   */
  vector<Lorentz5Momentum> _qnewlab;

  /**
   *  Momenta of the photons in the lab frame
   */
  vector<Lorentz5Momentum> _llab;

  /**
   * Total momentum of the photons in the lab frame
   */
  Lorentz5Momentum _bigLlab;
  //@}


  /**
   *  Reweighting factors due to differences between the true and crude
   *  distributions
   */
  //@{
  /**
   *  Reweighting factor for the real emission
   */
  double _dipolewgt;

  /**
   *  Reweighting factor for the YFS form-factor
   */
  double _yfswgt;

  /**
   *  Reweighting factor due to phase space
   */
  double _jacobianwgt;

  /**
   *  Reweighting factor due to matrix element corrections
   */
  double _mewgt;

  /**
   *  Maximum weight
   */
  double _maxwgt;
  //@}

  /**
   *  Angles of the photons with respect to the first charged particle
   * which are stored for numerical accuracy
   */
  //@{
  /**
   *  Cosine of the photon angles
   */
  vector<double> _cosphot;

  /**
   *  Sine of the photon angles
   */
  vector<double> _sinphot;
  //@}

  /**
   *  Weights for the individual photons
   */
  vector<double> _photonwgt;

  /**
   *  Whether a given photon passes the energy cut
   */
  vector<bool> _photcut;

  /**
   *  Type of unweighting to perform
   */
  unsigned int _mode;

  /**
   *  Maximum number of attempts to generate a result
   */
  unsigned int _maxtry;

  /**
   *  Option for the energy cut-off
   */
  unsigned int _energyopt;

  /**
   *  Option for the inclusion of higher order corrections
   */
  unsigned int _betaopt;

  /**
   *  Option for the form of the primary distribution
   */
  unsigned int _dipoleopt;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FFDipole. */
template <>
struct BaseClassTrait<Herwig::FFDipole,1> {
  /** Typedef of the first base class of FFDipole. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FFDipole class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FFDipole>
  : public ClassTraitsBase<Herwig::FFDipole> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::FFDipole"; }
};

/** @endcond */

}

#include "FFDipole.icc"

#endif /* HERWIG_FFDipole_H */
