// -*- C++ -*-
//
// SemiLeptonicScalarDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SemiLeptonicScalarDecayer_H
#define HERWIG_SemiLeptonicScalarDecayer_H
//
// This is the declaration of the SemiLeptonicScalarDecayer class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/FormFactors/ScalarFormFactor.h"
#include "Herwig/Decay/WeakCurrents/LeptonNeutrinoCurrent.h"
#include "ThePEG/Helicity/LorentzTensor.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>SemiLeptonicScalarDecayer</code> class is designed for the
 *  semileptonic decay of a (pesudo)scalar meson to another meson and a
 *  lepton/neutino pair.
 *
 *  This class implements the matrix element in a general form with form-factors
 *  for the hadronic current. These form-factor are specified using a ScalarFormFactor
 *  class and the leptonic part of the decay uses the LeptonNeutrinoCurrent.
 * 
 * @see DecayIntegrator
 * @see ScalarFormFactor
 * @see LeptonNeutrinoCurrent 
 *
 */
class SemiLeptonicScalarDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  SemiLeptonicScalarDecayer();
  
  /**
   * Check if this decayer can perfom the decay for a particular mode.
   * Uses the modeNumber member but can be overridden
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const ;

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
  static ClassDescription<SemiLeptonicScalarDecayer> initSemiLeptonicScalarDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  SemiLeptonicScalarDecayer & operator=(const SemiLeptonicScalarDecayer &);

private:

  /**
   * The leptonic current
   */
  LeptonNeutrinoCurrentPtr _current;

  /**
   * The form factor
   */
  ScalarFormFactorPtr _form;

  /**
   * the maximum weight for the integration of a given decay
   */
  vector<double> _maxwgt;

  /**
   * mapping of the mode to the form factor
   */
  vector<int> _modemap;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization vectors for the decay products
   */
  mutable vector<Helicity::LorentzPolarizationVector> _vectors;

  /**
   *  Polarization vectors for the decay products
   */
  mutable vector<Helicity::LorentzTensor<double> > _tensors;

  /**
   *  Constants for the mapping of the leptonic current
   */
  mutable vector<unsigned int> _constants;

  /**
   *  Spins of the particles
   */
  mutable vector<PDT::Spin> _ispin;

  /**
   *  Location of the outgoing meson
   */
  mutable unsigned int _imes;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of SemiLeptonicScalarDecayer.
 */
template <>
struct BaseClassTrait<Herwig::SemiLeptonicScalarDecayer,1> {
    /** Typedef of the base class of SemiLeptonicScalarDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::SemiLeptonicScalarDecayer>
  : public ClassTraitsBase<Herwig::SemiLeptonicScalarDecayer> {
  /** Return the class name. */
  static string className() { return "Herwig::SemiLeptonicScalarDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwSMDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_SemiLeptonicScalarDecayer_H */
