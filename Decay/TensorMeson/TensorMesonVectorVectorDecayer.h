// -*- C++ -*-
//
// TensorMesonVectorVectorDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TensorMesonVectorVectorDecayer_H
#define HERWIG_TensorMesonVectorVectorDecayer_H
//
// This is the declaration of the TensorMesonVectorVectorDecayer class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzTensor.h"

namespace Herwig {
using namespace ThePEG; 

/** \ingroup Decay
 *
 *  The TensorMesonVectorVectorDecayer class is designed to simulate
 *  the decay of a tensor meson to two spin-1 particles. In practice, at least
 *  for the light tensor mesons, this is only the decay of a tensor meson to two
 *  photons. In principle for bottom and charm tensors this may be the decay to
 *  two vector mesons an d it is used to simulate three and four body 
 *  \f$a_2\f$, \f$f_2\f$ and \f$K_2\f$ decays. 
 *  
 * The matrix element is
 *  given by \f[\mathcal{M}=\epsilon_0^{\alpha_1\alpha_2} \epsilon_1^{\beta} \epsilon_2^{\gamma} 
 *  \left(g_{\alpha_1\beta}+\frac{p_{1\alpha_1}p_{0\beta}}{p_0\cdot p_1-m_0m_1}\right)
 *  \left(g_{\alpha_2\gamma}+\frac{p_{2\alpha_2}p_{0\gamma}}{p_0\cdot p_2-m_0m_2}\right)\f]
 *  where \f$p_{0,1,2}\f$ are the momenta of the incoming pseudotensor and outgoing vector mesons, respectively,
 * and \f$\epsilon_{0}\f$ is the polarization tensor of the incoming pseudotensor meson and $\epslion_{1,2}$ are those of the outgoing
 * vector mesons.
 *
 *  The incoming tensor mesons together with their decay products and the coupling 
 *  \f$g\f$ can be specified using the interfaces for the class. The maximum weights
 *  for the decays can be calculated using the Initialize interface of the
 *  DecayIntegrator class or specified using the interface.
 *
 * @see DecayIntegrator
 * * @see \ref TensorMesonVectorVectorDecayerInterfaces "The interfaces"
 * defined for TensorMesonVectorVectorDecayer.
 *
 */
class TensorMesonVectorVectorDecayer: public DecayIntegrator {

public:

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
   *               in the GenericWidthGenerator class, in this case 9.
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
  TensorMesonVectorVectorDecayer & operator=(const TensorMesonVectorVectorDecayer &) = delete;

public:

  /**
   *   Set the parameters for a decay mode
   */
  string setUpDecayMode(string arg);

private:

  /**
   * the PDG codes for the incoming particles
   */
  vector<int> incoming_;

  /**
   * the PDG codes for the outgoing particles
   */
  vector<pair<int,int> > outgoing_;

  /**
   * the coupling for the decay
   */
  vector<Energy> coupling_;

  /**
   * the maximum weight for the decay
   */
  vector<double> maxWeight_;

  /**
   *  Storage of polarization tensors to try and increase
   *  speed
   */
  mutable vector<Helicity::LorentzTensor<double> > tensors_;

  /**
   *  Storage of the polarization vectors 
   */
  mutable vector<Helicity::LorentzPolarizationVector > vectors_[2];

  /**
   *   Storage of the \f$\rho\f$ matrix
   */
  mutable RhoDMatrix rho_;

};

}

#endif /* HERWIG_TensorMesonVectorVectorDecayer_H */
