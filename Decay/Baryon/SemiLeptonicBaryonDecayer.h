// -*- C++ -*-
#ifndef HERWIG_SemiLeptonicBaryonDecayer_H
#define HERWIG_SemiLeptonicBaryonDecayer_H
//
// This is the declaration of the SemiLeptonicBaryonDecayer class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/WeakCurrents/LeptonNeutrinoCurrent.h"
#include "Herwig/Decay/FormFactors/BaryonFormFactor.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzRSSpinorBar.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>SemiLeptonicBaryonDecayer</code> class is designed for the
 *  semi-leptonic decay of the baryons. It combine the form-factors from a 
 *  class inheriting from the BaryonFormFactor class and the leptonic current.
 *
 *  The decays of spin-\f$\frac12\f$ baryons to spin-\f$\frac12\f$ and
 *  spin-\f$\frac32\f$ baryons are currently supported. The only other decays
 *  which seem to occur in nature is the semi-leptonic decay of the \f$\Omega^-\f$
 *  which is \f$\frac32\to\frac12\f$.
 *
 * @see BaryonFormFactor.
 * 
 */
class SemiLeptonicBaryonDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  SemiLeptonicBaryonDecayer();

  /**
   * Check if this decayer can perfom the decay for a particular mode.
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
			 const tPDVector & children) const;
  
  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;
  
  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt The option for the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const ParticleVector & decay,MEOption meopt) const;

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

protected:

  /**
   * Matrix element for \f$\frac12\to\frac12\f$.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param inpart The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt The option for the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double halfHalf(const int ichan,const Particle & inpart,
		  const ParticleVector & decay,MEOption meopt) const;

  /**
   * Matrix element for \f$\frac12\to\frac32\f$.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param inpart The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt The option for the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double halfThreeHalf(const int ichan,const Particle & inpart,
		       const ParticleVector & decay,MEOption meopt) const;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<SemiLeptonicBaryonDecayer> initSemiLeptonicBaryonDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  SemiLeptonicBaryonDecayer & operator=(const SemiLeptonicBaryonDecayer &);

private:

  /**
   * The current for the leptons
   */
  LeptonNeutrinoCurrentPtr _current;

  /**
   * form-factor
   */
  BaryonFormFactorPtr _form;

  /**
   * the maximum weight
   */
  vector<double> _maxwgt;

  /**
   * mapping of the mode to the form-factor
   */
  vector<int> _modemap; 

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *   Spin-\f$\frac12\f$ spinors
   */
  mutable vector<LorentzSpinor<SqrtEnergy> > _inHalf;

  /**
   *   Spin-\f$\frac12\f$ barred spinors
   */
  mutable vector<LorentzSpinorBar<SqrtEnergy> > _inHalfBar;

  /**
   *   Spin-\f$\frac32\f$ spinors
   */
  mutable vector<LorentzRSSpinor<SqrtEnergy> > _inThreeHalf;

  /**
   *   Spin-\f$\frac32\f$ barred spinors
   */
  mutable vector<LorentzRSSpinorBar<SqrtEnergy> > _inThreeHalfBar;

  /**
   *  Constants for the mapping of the leptonic vector 
   */
  mutable vector<unsigned int> _constants;

  /**
   *  Spins of the particles
   */
  mutable vector<PDT::Spin> _ispin;

  /**
   *  Location of the outgoing baryon
   */
  mutable unsigned int _ibar;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of SemiLeptonicBaryonDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::SemiLeptonicBaryonDecayer,1> {
   /** Typedef of the base class of  SemiLeptonicBaryonDecayer. */
   typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::SemiLeptonicBaryonDecayer>
  : public ClassTraitsBase<Herwig::SemiLeptonicBaryonDecayer> {
   /** Return the class name.*/
  static string className() { return "Herwig::SemiLeptonicBaryonDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwBaryonDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_SemiLeptonicBaryonDecayer_H */
