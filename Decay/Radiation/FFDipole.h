// -*- C++ -*-
//
// FFDipole.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_FFDipole_H
#define HERWIG_FFDipole_H
//
// This is the declaration of the FFDipole class.
//

#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Decay/DecayIntegrator.fh"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Utilities/Maths.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/Interface/Interfaced.h"
#include "FFDipole.fh"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Constants::pi;

/** \ingroup Decay
 *
 * The FFDipole class generates radiation from a final-final dipole for
 * the generation of photons in decay by the SOPTHY algorithm.
 * 
 * @see SOPTHY 
 * @see \ref FFDipoleInterfaces "The interfaces"
 * defined for FFDipole.

 */
class FFDipole: public Interfaced {

public:

  /**
   * The default constructor.
   */
  FFDipole() :
    _emin(1.e-6*MeV), _eminrest(100*MeV), _eminlab(100*MeV), _emax(),
    _multiplicity(), _m(3), _charge(), _qdrf(2),
    _qnewdrf(2), _qprf(2), _qnewprf(2), _qlab(2), _qnewlab(2), _dipolewgt(),
    _yfswgt(), _jacobianwgt(), _mewgt(), _maxwgt(7.0), _mode(1), _maxtry(500),
    _energyopt(1), _betaopt(4), _dipoleopt(), _nweight(0), _wgtsum(0.), _wgtsq(0.),
    _weightOutput(false) {}

  /**
   *  Destructor
   */
  virtual ~FFDipole();

public:

  /**
   * Member to generate the photons from the dipole
   * @param p The decaying particle
   * @param children The decay products
   * @param decayer The decayer for this mode
   * @return The decay products with additional radiation
   */
  virtual ParticleVector generatePhotons(const Particle & p,
					 ParticleVector children,
					 tDecayIntegratorPtr decayer);

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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

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
  double exactDipoleWeight(double beta1,double ombeta1,
			   double beta2,double ombeta2,unsigned int iphot) {
    double opbc,ombc;
    // if cos is greater than zero use result accurate as cos->1
    if(_cosphot[iphot]>0) {
      opbc=1.+beta2*_cosphot[iphot];
      ombc=ombeta1+beta1*sqr(_sinphot[iphot])/(1.+_cosphot[iphot]);
    }
    // if cos is less    than zero use result accurate as cos->-1
    else {
      opbc=ombeta2+beta2*sqr(_sinphot[iphot])/(1.-_cosphot[iphot]);
      ombc=1.-beta1*_cosphot[iphot];
    }
    return 0.5/(opbc*ombc)*(1.+beta1*beta2
			    -0.5*ombeta1*(1.+beta1)*opbc/ombc		 
			    -0.5*ombeta2*(1.+beta2)*ombc/opbc);
  }
  
  /**
   * Jacobian factor for the weight
   */
  double jacobianWeight() {
    Energy pcm1=Kinematics::pstarTwoBodyDecay(_m[0],_m[1],_m[2]);
    Energy m12 =sqrt((_qnewdrf[0]+_qnewdrf[1]).m2())            ;
    Energy pcm2=Kinematics::pstarTwoBodyDecay(m12,_m[1],_m[2])  ;
    double betaprobeta = pcm2*_m[0]/pcm1/m12   ;
    double spros       = sqr(m12/_m[0])        ;
    double deltafn     = m12/(m12+_bigLdrf.e());
    return betaprobeta*spros*deltafn           ;
  }

  /**
   * Matrix element weight
   */
  double meWeight(const ParticleVector & children);

  /**
   * Member which generates the photons
   * @param boost Boost vector to take the particles produced back from
   * the decaying particle's rest frame to the lab
   * @param children The decay products
   */
  double makePhotons(const Boost & boost, 
		     const ParticleVector & children);

  /**
   *  Boost all the momenta from the dipole rest frame via the parent rest frame
   * to the lab
   * @param boost The boost vector from the rest frame to the lab
   * @return Whether or not it suceeded
   */
  bool boostMomenta(const Boost & boost);

  /**
   *  Remove any photons which fail the energy cuts
   * @return Number of photons removed
   */
  unsigned int removePhotons();

  /**
   *  The real emission weight in the collinear limit
   */
  double collinearWeight(const ParticleVector & children);

  /**
   *  The vrtiual correction weight
   */
  double virtualWeight(const ParticleVector & children);

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

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
  FFDipole & operator=(const FFDipole &) = delete;

private:

  /**
   * Debug output
   **/
  void printDebugInfo(const Particle & p,
		      const ParticleVector & children,
		      double wgt) const;

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

  /**
   *  The decayer
   */
  tDecayIntegratorPtr _decayer;

  /**
   *  The decaying particle
   */
  tPPtr _parent;

  /**
   *  Storage of averages etc for testing
   */
  //@{
  /**
   *  Number of attempts
   */
  unsigned int _nweight;

  /**
   *  Sum of weights
   */
  double _wgtsum;

  /**
   *  Sum of squares of weights
   */
  double _wgtsq;

  /**
   *  Whether or not to output the averages
   */
  bool _weightOutput;
  //@}
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
  static string className() { return "Herwig::FFDipole"; }
  /** Return the name of the shared library be loaded to get
   *  access to the DecayRadiationGenerator class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwSOPHTY.so"; }
};

/** @endcond */

}

#endif /* HERWIG_FFDipole_H */
