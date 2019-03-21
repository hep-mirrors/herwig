// -*- C++ -*-
//
// VectorMesonPVectorPScalarDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VectorMesonPVectorPScalarDecayer_H
#define HERWIG_VectorMesonPVectorPScalarDecayer_H
// This is the declaration of the VectorMesonPVectorPScalarDecayer class.

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  This class is designed for the decay of a vector meson to a pseudovector
 *  meson and a pseudoscalar meson.
 *  The current for the decay is
 *
 *  \f[\mathcal{M}= g\left[ p_V \cdot p_0 \epsilon_V\cdot \epsilon_0  
 *                     -p_V\cdot \epsilon_0 \epsilon_V \cdot p_0\right]\f]
 *
 * @see DecayIntegrator
 * @see \ref VectorMesonPVectorPScalarDecayerInterfaces "The interfaces"
 * defined for VectorMesonPVectorPScalarDecayer.
 * 
 *  \author Peter Richardson
 */
class VectorMesonPVectorPScalarDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  VectorMesonPVectorPScalarDecayer();

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

  /**
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class, in this case 4.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  bool twoBodyMEcode(const DecayMode & dm, int & mecode, double & coupling) const;

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
  static ClassDescription<VectorMesonPVectorPScalarDecayer> initVectorMesonPVectorPScalarDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  VectorMesonPVectorPScalarDecayer & operator=(const VectorMesonPVectorPScalarDecayer &) = delete;

private:

  /**
   * coupling for a decay
   */
  vector<InvEnergy> _coupling;

  /**
   * PDG codes for the incoming particles
   */
  vector<int> _incoming;

  /**
   * PDG codes for the outgoing pseudo-vector
   */
  vector<int> _outgoingA;

  /**
   * PDG codes for the outgoing pseudoscalar mesons.
   */
  vector<int> _outgoingP;

  /**
   * maximum weight for a decay
   */
  vector<double> _maxweight;

  /**
   *  Initial size of the vectors
   */
  unsigned int _initsize;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization vectors for the incoming and outgoing ferimons
   */
  mutable vector<Helicity::LorentzPolarizationVector> _vectors[2];

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of VectorMesonPVectorPScalarDecayer.
 */
template <>
struct BaseClassTrait<Herwig::VectorMesonPVectorPScalarDecayer,1> {
    /** Typedef of the base class of VectorMesonPVectorPScalarDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::VectorMesonPVectorPScalarDecayer>
  : public ClassTraitsBase<Herwig::VectorMesonPVectorPScalarDecayer> {
  /** Return the class name. */
  static string className() { return "Herwig::VectorMesonPVectorPScalarDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwVMDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_VectorMesonPVectorPScalarDecayer_H */
