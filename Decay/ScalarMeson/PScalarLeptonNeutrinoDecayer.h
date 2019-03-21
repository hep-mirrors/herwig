// -*- C++ -*-
//
// PScalarLeptonNeutrinoDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PScalarLeptonNeutrinoDecayer_H
#define HERWIG_PScalarLeptonNeutrinoDecayer_H
// This is the declaration of the PScalarLeptonNeutrinoDecayer class.

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzSpinorBar.h"

namespace Herwig {
using namespace ThePEG;

/**  \ingroup Decay
 *
 * The PScalarLeptonNeutrinoDecayer class is designed for the decay of 
 * pseudoscalar mesons to a lepton and a neutrino. Although it can be used
 * for charged pion and kaon decays it is mainly intended for the leptonic
 * decays of bottom and charm mesons.
 *
 *  The matrix element is given by 
 * \f[\mathcal{M} = \frac1{\sqrt{2}}f_PG_FV_{CKM}m_l\bar{u}(p_{\ell})(1-\gamma_5)v(p_\nu),\f]
 * where
 * - \f$f_P\f$ is the pseudoscalar decay constant.
 * - \f$G_F\f$ is the Fermi constant
 * - \f$V_{CKM}\f$ is the relevant CKM matrix element
 * - \f$p_\ell\f$ is the momentum of the charged lepton
 * - \f$p_\nu\f$ is the momentum of the neutrino
 *
 * @see DecayIntegrator
 * 
 */
class PScalarLeptonNeutrinoDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  PScalarLeptonNeutrinoDecayer();
  
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
  double me2( const int ichan,const Particle & part,
	     const ParticleVector & decay, MEOption meopt) const;

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
  static ClassDescription<PScalarLeptonNeutrinoDecayer> initPScalarLeptonNeutrinoDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  PScalarLeptonNeutrinoDecayer & operator=(const PScalarLeptonNeutrinoDecayer &) = delete;

private:

  /**
   * the PDG code for the incoming particle
   */
  vector<int> _incoming;

  /**
   * the meson decay constant for a particular particle multiplied by the CKM matrix
   * element, \e i.e. \f$f_pV_{CKM}\f$
   */
  vector<Energy> _decayconstant;

  /**
   * which outgoing leptons are allowed for a particular decay
   */
  vector<unsigned int> _leptons;

  /**
   * the maximum weight for the integration of a given decay to \f$e\nu_e\f$.
   */
  vector<double> _maxweighte;

  /**
   * the maximum weight for the integration of a given decay to \f$\mu\nu_\mu\f$.
   */
  vector<double> _maxweightmu;

  /**
   * the maximum weight for the integration of a given decay to \f$\tau\nu_\tau\f$.
   */
  vector<double> _maxweighttau;

  /**
   *  initial number of modes
   */
  unsigned int _initsize;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Spinors for the decay products
   */
  mutable vector<Helicity::LorentzSpinor   <SqrtEnergy> > _wave;

  /**
   *  barred spinors for the decay products
   */
  mutable vector<Helicity::LorentzSpinorBar<SqrtEnergy> > _wavebar;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

template <>
/**
 * The following template specialization informs ThePEG about the
 * base class of PScalarLeptonNeutrinoDecayer.
 */
struct BaseClassTrait<Herwig::PScalarLeptonNeutrinoDecayer,1> {
    /** Typedef of the base class of PScalarLeptonNeutrinoDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::PScalarLeptonNeutrinoDecayer>
  : public ClassTraitsBase<Herwig::PScalarLeptonNeutrinoDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig::PScalarLeptonNeutrinoDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwSMDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_PScalarLeptonNeutrinoDecayer_H */
