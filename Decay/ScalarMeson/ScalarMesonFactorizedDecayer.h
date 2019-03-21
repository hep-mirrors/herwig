// -*- C++ -*-
//
// ScalarMesonFactorizedDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ScalarMesonFactorizedDecayer_H
#define HERWIG_ScalarMesonFactorizedDecayer_H
//
// This is the declaration of the ScalarMesonFactorizedDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/WeakCurrents/WeakDecayCurrent.h"
#include "Herwig/Decay/FormFactors/ScalarFormFactor.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>ScalarMesonFactorizedDecayer</code> class is a class which combines a 
 * WeakDecayCurrent and a ScalarFormFactor in the naive factorization approximation
 * to perform the non-leptonic weak decays of scalar mesons.
 *
 * @see DecayIntegrator
 * @see WeakDecayCurrent
 * @see ScalarFormFactor
 *
 */

class ScalarMesonFactorizedDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  ScalarMesonFactorizedDecayer();

public:

  /** @name Virtual functions required by the Decayer and DecayIntegrator classes. */
  //@{
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;
  
  /**
   * Check if this decayer can perfom the decay for a particular mode.
   * Uses the modeNumber member but can be overridden
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * This function combines the current and the form factor to give the matrix
   * element.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2( const int ichan, const Particle & part,
		     const ParticleVector & decay, MEOption meopt) const;
  //@}

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
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans)
   ;

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}


private:

  /**
   * Find duplicate modes in the list of particles
   * @param imode The mode we are studying
   * @param particles The external particles for the different modes
   * @param loc The location of the duplicate mode
   * @param cc  If the duplicate is the charge conjugate
   */
  void findModes(unsigned int imode,vector<tPDVector> & particles,
		 vector<unsigned int> & loc,vector<bool> & cc);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ScalarMesonFactorizedDecayer> initScalarMesonFactorizedDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ScalarMesonFactorizedDecayer & operator=(const ScalarMesonFactorizedDecayer &) = delete;

private:

  /**
   *  The weak decay current
   */
  vector<WeakDecayCurrentPtr> _current;

  /**
   * The baryon form factor
   */
  vector<ScalarFormFactorPtr> _form;

  /**
   *  The perturbative coefficients
   */
  //@{
  /**
   *  The perturbative \f$a_1\f$ coefficient for b decays.
   */
  double _a1b;

  /**
   *  The perturbative \f$a_2\f$ coefficient for b decays.
   */
  double _a2b;

  /**
   *  The perturbative \f$a_1\f$ coefficient for c decays.
   */
  double _a1c;

  /**
   *  The perturbative \f$a_2\f$ coefficient for c decays.
   */
  double _a2c;
  //@}

  /**
   * Mapping of the modes to the currents
   */
  //@{
  /**
   *  First map
   */
  vector<vector<unsigned int> > _currentmapA;

  /**
   *  Second map
   */
  vector<vector<unsigned int> > _currentmapB;
  //@}

  /**
   * Mapping of the modes to the form factors
   */
  //@{
  /**
   *  First map
   */
  vector<vector<unsigned int> > _formmapA;

  /**
   *  Second map
   */
  vector<vector<unsigned int> > _formmapB;
  //@}
  /**
   *  Outgoing particle from the form factor
   */
  vector<vector<unsigned int> > _formpart;

  /**
   * The CKM factors
   */
  vector<vector<Complex> > _CKMfact;

  /**
   * location of the weights
   */
  vector<int> _wgtloc;

  /**
   * the maximum weights
   */
  vector<double> _wgtmax;

  /**
   *  Weights for the different channels
   */
  vector<double> _weights;

  /**
   * Pointer to the CKM object.
   */
  Ptr<StandardCKM>::pointer _ckm;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization vectors for the decay products
   */
  mutable vector<vector<Helicity::LorentzPolarizationVector> > _vectors;

  /**
   *  Polarization tensors for the decay products
   */
  mutable vector<vector<Helicity::LorentzTensor<double>    > > _tensors;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ScalarMesonFactorizedDecayer. */
template <>
 struct BaseClassTrait<Herwig::ScalarMesonFactorizedDecayer,1> {
  /** Typedef of the first base class of ScalarMesonFactorizedDecayer. */
   typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ScalarMesonFactorizedDecayer class and the shared object where it is defined. */
template <>
 struct ClassTraits<Herwig::ScalarMesonFactorizedDecayer>
  : public ClassTraitsBase<Herwig::ScalarMesonFactorizedDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ScalarMesonFactorizedDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the ScalarMesonFactorizedDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwSMDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ScalarMesonFactorizedDecayer_H */
