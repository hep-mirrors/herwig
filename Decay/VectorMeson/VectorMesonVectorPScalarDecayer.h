// -*- C++ -*-
//
// VectorMesonVectorPScalarDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VectorMesonVectorPScalarDecayer_H
#define HERWIG_VectorMesonVectorPScalarDecayer_H
// This is the declaration of the VectorMesonVectorPScalarDecayer class.

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {

/** \ingroup Decay
 *
 *  This class is designed for the decay of a vector meson to another spin-1 
 *  particle, either another vector meson or a photon, and a pseduoscalar meson.
 *  The matrix element for the decay is 
 *  \f[\mathcal{M} = g\epsilon^{\mu\nu\alpha\beta} 
 *     \epsilon_{0\mu} p_{0\nu}  p_{1\alpha} \epsilon_{1\beta}, \f]
 *  where \f$p_0\f$ is the momentum of the decaying particle, \f$p_1\f$ is the momentum
 *  of the outgoing vector particle, \f$\epsilon_0\f$ is the polarization vector of
 * the incoming meson and \f$\epsilon_1\f$ is the polarization vector of the
 * outgoing vector particle..
 *
 *  Examples of such decays are \f$\rho\to\pi\gamma\f$.
 *
 *  The incoming vector mesons together with their decay products and the coupling 
 *  \f$g\f$ can be specified using the interfaces for the class. The maximum weights
 *  for the decays can be calculated using the Initialize interface of the
 *  DecayIntegrator class or specified using the interface.
 *
 *  The incoming and outgoing particles, couplings and maximum weights for
 *  many of the common \f$V\to PV\f$ decays are specified in the default
 *  constructor.
 *
 * @see DecayIntegrator
 * @see \ref VectorMesonVectorPScalarDecayerInterfaces "The interfaces"
 * defined for VectorMesonVectorPScalarDecayer.
 * 
 *  \author Peter Richardson
 *
 */
class VectorMesonVectorPScalarDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  VectorMesonVectorPScalarDecayer();

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
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const tPDVector & outgoing,
	     const vector<Lorentz5Momentum> & momenta,
	     MEOption meopt) const;

  /**
   *   Construct the SpinInfos for the particles produced in the decay
   */
  virtual void constructSpinInfo(const Particle & part,
				 ParticleVector outgoing) const;

  /**
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class, in this case 1.
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
   * Private and non-existent assignment operator.
   */
  VectorMesonVectorPScalarDecayer & operator=(const VectorMesonVectorPScalarDecayer &) = delete;

public:

  /**
   *   Set the parameters for a decay mode
   */
  string setUpDecayMode(string arg);

private:

  /**
   * coupling for a decay
   */
  vector<InvEnergy> coupling_;

  /**
   * the PDG codes for the incoming particle
   */
  vector<long> incoming_;

  /**
   * the PDG codes for the outgoing particles (vector first)
   */
  vector<pair<long,long> > outgoing_;

  /**
   * maximum weight for a decay
   */
  vector<double> maxweight_;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Polarization vectors
   */
  mutable vector<Helicity::LorentzPolarizationVector> vectors_[2];

};

}


#endif /* HERWIG_VectorMesonVectorPScalarDecayer_H */
