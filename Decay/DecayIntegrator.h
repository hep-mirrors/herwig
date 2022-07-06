// -*- C++ -*-
//
// DecayIntegrator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_DecayIntegrator_H
#define Herwig_DecayIntegrator_H
//
// This is the declaration of the DecayIntegrator class.
//

#include "DecayIntegrator.fh"
#include "HwDecayerBase.h"
#include "PhaseSpaceMode.fh"
#include "Herwig/PDT/WidthCalculatorBase.fh"
#include "Radiation/DecayRadiationGenerator.h"
#include <Herwig/Decay/DecayVertex.h>

namespace Herwig {

using namespace ThePEG;

  /** \ingroup Decay
   * \class DecayIntegrator
   * \brief Main class for Decayers implementing multi-channel phase space integration.
   * \author Peter Richardson
   *
   *  This class is designed to be the base class for Herwig decays including
   *  the implementation of a multichannel decayer or n-body phase space decays.
   *
   *  The <code>DecayIntegrator</code> class inherits from ThePEG's Decayer class
   *  and makes use of the <code>PhaseSpaceMode</code> class to specify a number
   *  of decay modes.
   *
   *  Additional modes can be added using the addMode method. In practice the 
   *  phase space channels for a particular mode are usually constructed in the 
   *  doinit member of a Decayer and then the modes added to the Decayer.
   *
   *  For the majority of the decays currently implemented the 
   *  phase-space integration has been optimised and the maximum weight set.
   *  If the parameters of the decay model are changed the Initialize interface 
   *  can be used to optimise the integration and calculate the maximum weight.
   *
   *  In classes inheriting from this the me2() member which gives the matrix element
   *  squared must be implemented. This should be combined with the setting of the
   *  phase space channels, and the setting of which channels to use and their
   *  initial weights in the doinit() member. The different decay modes should then
   *  be initialized in the initrun() member if needed. The generate member can then
   *  be called from the decay() member to generate a phase-space configuration for a 
   *  decay.
   *   
   * @see DecayPhaseSpaceMode
   * @see DecayPhaseSpaceChannel
   * @see \ref DecayIntegratorInterfaces "The interfaces"
   * defined for DecayIntegrator.
   */
class DecayIntegrator: public HwDecayerBase {

public:
  
  /**
   *  and DecayPhaseMode
   */
  friend class PhaseSpaceMode;

  /**
   *  Enum for the matrix element option
   */
  enum MEOption {Initialize,Calculate,Terminate};

public:

  /**
   * The default constructor.
   */
  DecayIntegrator() : nIter_(10), nPoint_(10000), nTry_(500),
		       generateInter_(false), iMode_(-1),
		       realME_(false), virtualME_(false), eps_(ZERO), warnings_(true)
  {}

public:
  
  /**
   * Check if this decayer can perfom the decay for a particular mode.
   * Uses the modeNumber member but can be overridden
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const {
    bool cc;
    return modeNumber(cc,parent,children)>=0;
  }
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const tPDVector & children) const;
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const = 0;

  /**
   * The mode being used for this decay
   */
  int imode() const {return iMode_;}

  /**
   * Add a phase-space mode to the list
   * @param mode The mode being added.
   */
  void addMode(PhaseSpaceModePtr mode) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int ichan, const Particle & part,
		     const tPDVector & outgoing,
		     const vector<Lorentz5Momentum> & momenta,
		     MEOption meopt) const = 0;

  /**
   *   Construct the SpinInfos for the particles produced in the decay
   */
  virtual void constructSpinInfo(const Particle & part,
				 ParticleVector outgoing) const = 0;
  
  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;
  
  /**
   * Set the code for the partial width. Finds the partial width in the
   * GenericWidthGenerator class which corresponds to the decay mode.
   * @param dm The DecayMode
   * @param imode The mode. 
   */
  void setPartialWidth(const DecayMode & dm, int imode);
  /**
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  virtual bool twoBodyMEcode(const DecayMode &, int & mecode,
			     double & coupling) const {
    coupling = 1.;
    mecode   = -1;
    return false;
  }				     
  
  /**
   * Method to return an object to calculate the 3 (or higher body) partial width
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;
  
  /**
   * The matrix element to be integrated for the three-body decays as a function
   * of the invariant masses of pairs of the outgoing particles.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element
   */
  virtual double threeBodyMatrixElement(const int imode,  const Energy2 q2,
					const Energy2 s3, const Energy2 s2, 
					const Energy2 s1, const Energy  m1, 
					const Energy  m2, const Energy  m3) const;
  
  /**
   * The differential three body decay rate with one integral performed.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s  The invariant mass which still needs to be integrate over.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The differential rate \f$\frac{d\Gamma}{ds}\f$
   */
  virtual InvEnergy threeBodydGammads(const int imode, const Energy2 q2,
				      const Energy2 s,
				      const Energy m1, const Energy m2, 
				      const Energy m3) const;

  /**
   * Finds the phase-space mode corresponding to a given decay mode
   * @param dm The DecayMode
   */
  int findMode(const DecayMode & dm);
  
public:

  /**
   *  Members for the generation of QED radiation in the decays
   */
  //@{
  /**
   * Use the DecayRadiationGenerator to generate photons in the decay.
   * @param p The Particle instance being decayed
   * @param children The decay products
   * @return A particle vector containing the decay products after the generation
   * of photons.
   */
  ParticleVector generatePhotons(const Particle & p,ParticleVector children) {
    return photonGen_->generatePhotons(p,children,this);
  }

  /**
   *  check if photons can be generated in the decay
   */
  bool canGeneratePhotons() {return photonGen_;}

  /**
   *  The one-loop virtual correction.
   * @param imode The mode required.
   * @param part  The decaying particle.
   * @param products The decay products including the radiated photon.
   * @return Whether the correction is implemented
   */
  virtual double oneLoopVirtualME(unsigned int imode,
				  const Particle & part, 
				  const ParticleVector & products);

  /**
   *  Whether or not the one loop matrix element is implemented
   */
  bool hasOneLoopME() {return virtualME_;}
  
  /**
   *  The real emission matrix element
   * @param imode The mode required
   * @param part  The decaying particle
   * @param products The decay products including the radiated photon
   * @param iemitter The particle which emitted the photon
   * @param ctheta   The cosine of the polar angle between the photon and the
   *                 emitter
   * @param stheta   The sine of the polar angle between the photon and the
   *                 emitter 
   * @param rot1 Rotation from rest frame to frame for real emission
   * @param rot2 Rotation to place emitting particle along z
   */
  virtual InvEnergy2 realEmissionME(unsigned int imode,
				    const Particle & part, 
				    ParticleVector & products,
				    unsigned int iemitter,
				    double ctheta, double stheta,
				    const LorentzRotation & rot1,
				    const LorentzRotation & rot2);

  /**
   *  Whether or not the real emission matrix element is implemented
   */
  bool hasRealEmissionME() {return realME_;}
  //@}
  
public:

  /**
   *  The output operator is a friend, this is mainly for debugging
   */    
  friend ostream & operator<<(ostream & os, const DecayIntegrator & decay);

  
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

protected:

  /**
   * Generate the momenta for the decay
   * @param inter Generate the intermediates produced in the decay as well as the
   * final particles.
   * @param cc Is this the mode defined or its charge conjugate.
   * @param imode The mode being generated.
   * @param inpart The decaying particle.
   * @return The particles produced inthe decay.
   */
  ParticleVector generate(bool inter,bool cc, const unsigned int & imode,
			  const Particle & inpart) const;

  /**
   * Set the mode being use for this decay.
   */
  void imode(int in) { iMode_ = in;}
  
  /**
   * Set the helicity matrix element for the decay.
   */
  void ME(DecayMEPtr in) const { matrixElement_ = in;}
  
  /**
   * The helicity amplitude matrix element for spin correlations.
   */
  DecayMEPtr ME() const {return matrixElement_;}
   
  /**
   * Reset the properities of all intermediates.
   * @param part The intermediate particle being reset.
   * @param mass The mass of the particle.
   * @param width The width of the particle.
   */
  void resetIntermediate(tcPDPtr part, Energy mass, Energy width);

  /**
   * Initialize the phase-space mode
   * @param imode The mode
   * @param init Whether or not to perform the initialization
   */
  Energy initializePhaseSpaceMode(unsigned int imode,bool init, bool onShell=false) const;
  
protected:
  
  /**
   *    Methods to set variables in inheriting classes
   */
  //@{
  /**
   * Set whether or not the intermediates are included 
   */ 
  void generateIntermediates(bool in) {generateInter_=in;}

  /**
   * Set whether or not the intermediates are included 
   */ 
  bool generateIntermediates() const {return generateInter_;}
  
  /**
   *  Whether or not the one loop matrix element is implemented
   */
  void hasOneLoopME(bool in) {virtualME_=in;}

  /**
   *  Whether or not the real emission matrix element is implemented
   */
  void hasRealEmissionME(bool in) {realME_=in;}

  /**
   * Set the epsilon parameter
   */
  void epsilonPS(Energy in) {eps_=in;}

  /**
   *  Clear the models
   */
  void clearModes() {modes_.clear();}

protected:
  
  /**
   * Number of decay modes
   */
  unsigned int numberModes() const {return modes_.size();}

  /**
   * Pointer to a mode
   */
  tPhaseSpaceModePtr mode(unsigned int ix) {
    return modes_[ix];
  }
  /**
   * Pointer to a mode
   */
  tcPhaseSpaceModePtr mode(unsigned int ix) const {
    return modes_[ix];
  }
  
public:

bool warnings() const {
  return warnings_;
}

private:

  /**
   * Private and non-existent assignment operator.
   */
  DecayIntegrator & operator=(const DecayIntegrator &) = delete;

  /**
   *  Parameters for the integration
   */
  //@{
  /**
   * Number of iterations for th initialization.
   */
  unsigned int nIter_;

  /**
   * Number of points for initialisation
   */
  unsigned int nPoint_;

  /**
   * number of attempts to generate the decay
   */
  unsigned int nTry_;
  
  /**
   * List of the decay modes
   */
  mutable vector<PhaseSpaceModePtr> modes_;

  //@}
  
  /**
   *  Whether to include the intermediates whne outputing the results.
   */
  bool generateInter_;
  
  /**
   *  Pointer to the object generating the QED radiation in the decay
   */
  DecayRadiationGeneratorPtr photonGen_;

  /**
   * mode currently being generated  
   */
  mutable int iMode_;

  /**
   * The helicity matrix element for the current decay
   */
  mutable DecayMEPtr matrixElement_;

  /**
   *  Whether or not the real photon emission matrix element exists
   */
  bool realME_;

  /**
   *  Whether or not the one-loop matrix element exists
   */
  bool virtualME_;

  /**
   *   Epsilon parameter for phase-space integration
   */
  Energy eps_;

  /**
   * option for turinh on/off log warnings in Phase class
   */
  bool warnings_;

protected:

  /**
   * Exception for this class and those inheriting from it
   */
  class DecayIntegratorError: public Exception {};

};
  
/**
 * Output information on the DecayIntegrator for debugging purposes
 */
ostream & operator<<(ostream &, const DecayIntegrator &);

}

#endif /* Herwig_DecayIntegrator_H */
