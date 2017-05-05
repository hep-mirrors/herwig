// -*- C++ -*-
//
// PScalarVectorVectorDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PScalarVectorVectorDecayer_H
#define HERWIG_PScalarVectorVectorDecayer_H
// This is the declaration of the PScalarVectorVectorDecayer class.
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>PScalarVectorVectorDecayer</code> class is designed to perform the decay 
 * of a pseudoscalar meson to two spin-1 particles. The majority of these decays
 * are of a light pseudoscalar meson to \f$\gamma\gamma\f$ where
 * including the matrix-element is unnessecary. However there are a small number
 * of decays, 
 * \e e.g. \f$\eta'\to\omega\gamma\f$
 * where the use of this decayer is required to get the correct correlations.
 *
 *  The matrix element is taken to be 
 * \f[\mathcal{M} = g\epsilon^{\mu\nu\alpha\beta}
 *                   p_{1\mu}   \epsilon_{1\nu}
 *                   p_{2\alpha}\epsilon_{2\beta},
 *   \f]
 *  where \f$p_{1,2}\f$ and \f$\epsilon_{1,2}\f$ are the momenta and polarzation
 *  vectors of the outgoing vectors.
 *
 *  The incoming pseudoscalar meson, the outgoing vectors and the coupling can
 *  be specified using the relevant interfaces.
 *
 * @see DecayIntegrator
 * @see PScalarVectorFermionsDecayer
 * @see PScalar4FermionsDecayer
 *
 * \author Peter Richardson
 * 
 */
class PScalarVectorVectorDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  PScalarVectorVectorDecayer();
  
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
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class, in this case 3.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  bool twoBodyMEcode(const DecayMode & dm, int & mecode, double & coupling) const;

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
  static ClassDescription<PScalarVectorVectorDecayer> initPScalarVectorVectorDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  PScalarVectorVectorDecayer & operator=(const PScalarVectorVectorDecayer &);

private:

  /**
   * the PDG code for the incoming particle
   */
  vector<int> _incoming;

  /**
   * the PDG code for the first outgoing particle
   */
  vector<int> _outgoing1;

  /**
   * the PDG code for the second outgoing particle
   */
  vector<int> _outgoing2;

  /**
   * the coupling for the decay, \f$g\f$.
   */
  vector<InvEnergy> _coupling;

  /**
   * the maximum weight for the decay
   */
  vector<double> _maxweight;

  /**
   *  initial number of modes
   */
  unsigned int _initsize;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization vectors for the decay products
   */
  mutable vector<Helicity::LorentzPolarizationVector> _vectors[2];
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of PScalarVectorVectorDecayer.
 */
template <>
struct BaseClassTrait<Herwig::PScalarVectorVectorDecayer,1> {
    /** Typedef of the base class of PScalarVectorVectorDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::PScalarVectorVectorDecayer>
  : public ClassTraitsBase<Herwig::PScalarVectorVectorDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig::PScalarVectorVectorDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwSMDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_PScalarVectorVectorDecayer_H */
