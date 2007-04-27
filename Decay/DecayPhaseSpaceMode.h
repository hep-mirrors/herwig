// -*- C++ -*-
#ifndef HERWIG_DecayPhaseSpaceMode_H
#define HERWIG_DecayPhaseSpaceMode_H
//
// This is the declaration of the DecayPhaseSpaceMode class.
//
#include "ThePEG/Interface/Interfaced.h"
#include "DecayPhaseSpaceMode.fh"
#include "DecayPhaseSpaceChannel.h"
#include "Herwig++/PDT/GenericWidthGenerator.h"
#include "Herwig++/PDT/GenericMassGenerator.h"
#include <Herwig++/Helicity/Correlations/DecayVertex.h>
#include "DecayIntegrator.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>DecayPhaseSpaceMode</code> class is designed to store a group
 * of phase-space channels for use by the DecayIntegrator class to 
 * generate the phase-space for a given decay mode.
 *
 *  Additional phase-space channels can be added using the addChannel member.
 * 
 *  In practice the modes are usually constructed together with the a number of
 *  <code>DecayPhaseSpaceChannel</code> objects. In classes inheriting from the
 *  DecayIntegrator class.
 *
 * @see DecayIntegrator
 * @see DecayPhaseSpaceChannel
 *
 * @author  Peter Richardson
 * 
 */
class DecayPhaseSpaceMode: public Interfaced {

  /**
   * A friend operator to allow the mode to be outputted for debugging purposes.
   */
  friend ostream & operator<<(ostream &, const DecayPhaseSpaceMode &);

  /**
   * DecayIntegrator is a friend to avoid making many of the phase space
   * generation members public.
   */
  friend class DecayIntegrator;

  /**
   * DecayPhaseSpaceChannel is a friend to avoid making many of the phase space
   * generation members public
   */
  friend class DecayPhaseSpaceChannel;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline DecayPhaseSpaceMode();

  /**
   * Constructor with a pointer to a <code>DecayPhaseIntegrator</code> and a vector
   * of particle data objects for the external particles. This 
   * is the constructor which should normally be used in decayers.
   * @param in The particle data objects for the external particles
   * @param intin A pointer to the DecayIntegrator class using this mode.
   */
  inline DecayPhaseSpaceMode(PDVector in, tcDecayIntegratorPtr intin);
  //@}

  /**
   * Access to the external particles.
   * @param ix The external particle required.
   * @return A pointer to the ParticleData object.
   */
  inline tcPDPtr externalParticles(int ix) const;

  /**
   * Number of external particles.
   * @return The number of external particles.
   */
  inline unsigned int numberofParticles() const;

  /**
   * Number of channels
   * @return The number of channels.
   */
  inline unsigned int numberChannels() const;

  /**
   * Add a new channel. 
   * @param channel A pointer to the new DecayPhaseChannel
   */
  inline void addChannel(DecayPhaseSpaceChannelPtr channel);

  /**
   * Reset the properties of one of the intermediate particles. Only a specific channel
   * is reset.
   * @param ichan The channel to reset.
   * @param part The ParticleData object of the particle to reset
   * @param mass The mass of the intermediate.
   * @param width The width of gthe intermediate.
   */
  inline void resetIntermediate(int ichan, tcPDPtr part, Energy mass, Energy width);

  /**
   * Reset the properties of one of the intermediate particles. All the channels 
   * are reset.
   * @param part The ParticleData object of the particle to reset
   * @param mass The mass of the intermediate.
   * @param width The width of gthe intermediate.
   */
  inline void resetIntermediate(tcPDPtr part, Energy mass, Energy width);

  /**
   * Get the maximum weight for the decay.
   * @return The maximum weight.
   */
  inline double maxWeight() const;

  /**
   * Set the maximum weight for the decay.
   * @param wgt The maximum weight. 
   */
  inline void setMaxWeight(double wgt) const;

  /**
   * Get the weight for a channel. This is the weight for the multi-channel approach.
   * @param ichan The channel.
   * @return The weight for the channel.
   */
  inline double channelWeight(unsigned int ichan) const;

  /**
   * Set the weights for the different channels.
   * @param in The weights for the different channels in the multi-channel approach.
   */
  inline void setWeights(const vector<double> & in) const;

protected:

  /** @name Set-up, Initialization and Access Members */
  //@{
  /**
   * Initialise the phase space.
   * @param init Perform the initialization.
   */
  void initializePhaseSpace(bool init);

  /**
   * Set the integration parameters
   * @param iter The number of iterations to use for initialization.
   * @param points The number of points to use for each iteration during initialization.
   * @param ntry The number of tries to generate a decay.
   */
  inline void setIntegrate(int iter,int points,int ntry);

  /**
   * Set the partial width to use for normalization. This is the partial width
   * in the WidthGenerator object.
   * @param in The partial width to use.
   */
  void setPartialWidth(int in);

  //@}


  /** @name Phase-Space Generation Members */
  //@{

  /**
   * Access to the matrix element from the decayer.
   * @param bin Generate the vertex information for spin correlations
   * @param ichan The channel, this is to allow the matrix element to be used to
   *              select the intermediates
   * @param inpart The incoming particle.
   * @param outpart The outgoing particles.
   */
  inline double me2(bool bin,const int ichan ,const Particle &inpart,
		    const ParticleVector &outpart) const;

  /**
   * Generate the decay.
   * @param intermediates Whether or not to generate the intermediate particle
   *                      in the decay channel.
   * @param cc Whether we are generating the mode specified or the charge 
   *           conjugate mode.
   * @param inpart The incoming particle.
   * @return The outgoing particles.
   */
  ParticleVector generate(bool intermediates,bool cc,const Particle & inpart) const;


  /**
   * Select which channel to use to output the particles.
   * @param inpart  The incoming particles.
   * @param outpart The outgoing particles.
   */
  inline int selectChannel(const Particle & inpart, ParticleVector & outpart) const;

  /**
   * Return the weight for a given phase-space point.
   * @param vertex Produce the SpinInfo objects for the spin correlations.
   * @param cc Whether we are generating the mode specified or the charge 
   *           conjugate mode.
   * @param ichan The channel to generate the weight for.
   * @param in The incoming particle.
   * @param particles The outgoing particles.
   * @return The weight.
   */
  Energy weight(bool vertex, bool cc,int & ichan, const Particle & in,
		ParticleVector & particles) const;
    
  /**
   * Return the weight and momenta for a flat phase-space decay.
   * @param cc Whether we are generating the mode specified or the charge 
   *           conjugate mode.
   * @param inpart The incoming particle.
   * @param outpart The outgoing particles.
   * @return The weight.
   */
  Energy flatPhaseSpace(bool cc,const Particle & inpart, ParticleVector & outpart) const;
  
  /**
   * Generate a phase-space point using multichannel phase space.
   * @param cc Whether we are generating the mode specified or the charge 
   *           conjugate mode.
   * @param ichan The channel to generate the weight for.
   * @param in The incoming particle.
   * @param particles The outgoing particles.
   * @return The weight.
   */
  Energy channelPhaseSpace(bool cc,int & ichan, const Particle & in, 
			   ParticleVector & particles) const;

  /**
   * Construct the vertex for spin corrections
   * @param in The incoming particle.
   * @param out The outgoing particles.
   */
  void constructVertex(const Particle & in, const ParticleVector & out) const;

  /**
   * Generate the masses of the external particles.
   * @param inmass The mass of the decaying particle.
   * @param wgt The weight for the masses.
   * @return The masses.
   */
  vector<Energy> externalMasses(Energy inmass,double & wgt) const;

  //@}

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<DecayPhaseSpaceMode> initDecayPhaseSpaceMode;

  /**
   * Private and non-existent assignment operator.
   */
  DecayPhaseSpaceMode & operator=(const DecayPhaseSpaceMode &);

 private:

  /**
   * pointer to the decayer
   */
  tcDecayIntegratorPtr _integrator;

  /**
   * pointers to the phase-space channels
   */
  vector<DecayPhaseSpaceChannelPtr> _channels;

  /**
   * the weights for the different channels
   */
  mutable vector<double> _channelwgts;

  /**
   * the maximum weight for the decay
   */
  mutable double _MaxWeight;

  /**
   * Number of iterations for the initialization.
   */
  int _niter;

  /**
   * Number of weights for each iteration of the initialization.
   */
  int _npoint;

  /**
   * Number of attempts to generate the decay
   */
  int _ntry;

  /**
   * External particles
   */
  PDVector _extpart;

  /**
   * Which of the partial widths of the incoming particle to use
   */
  int _partial;

  /**
   * The width generator for the incoming particle.
   */
  cGenericWidthGeneratorPtr _widthgen;

  /**
   *  The mass generators for the outgoing particles.
   */
  vector<cGenericMassGeneratorPtr> _massgen;
};

  /**
   *  The output operator which is used for debugging.
   */
ostream & operator<<(ostream &, const DecayPhaseSpaceMode &);

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

template <>
/**
 * The following template specialization informs ThePEG about the
 * base class of DecayPhaseSpaceMode.
 */
 struct BaseClassTrait<Herwig::DecayPhaseSpaceMode,1> {
  /** Typedef of the base class of DecayPhaseSpaceMode. */
  typedef Interfaced NthBase;
};

template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
 struct ClassTraits<Herwig::DecayPhaseSpaceMode>
  : public ClassTraitsBase<Herwig::DecayPhaseSpaceMode> {
   /** Return the class name. */
   static string className() { return "Herwig++::DecayPhaseSpaceMode"; }
};

/** @endcond */

}
#include "DecayPhaseSpaceMode.icc"

#endif /* HERWIG_DecayPhaseSpaceMode_H */
