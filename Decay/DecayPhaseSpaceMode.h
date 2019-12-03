// -*- C++ -*-
//
// DecayPhaseSpaceMode.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DecayPhaseSpaceMode_H
#define HERWIG_DecayPhaseSpaceMode_H
//
// This is the declaration of the DecayPhaseSpaceMode class.
//
#include "ThePEG/Interface/Interfaced.h"
#include "DecayPhaseSpaceMode.fh"
#include "DecayPhaseSpaceChannel.h"
#include "Herwig/PDT/GenericWidthGenerator.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include <Herwig/Decay/DecayVertex.h>
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
class DecayPhaseSpaceMode: public Base {

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
  DecayPhaseSpaceMode() :  _maxweight(0.),_niter(10), _npoint(10000), _ntry(500),
			   _partial(-1), _testOnShell(false), _ichannel(999) {}

  /**
   * Constructor with a pointer to a <code>DecayPhaseIntegrator</code> and a vector
   * of particle data objects for the external particles. This 
   * is the constructor which should normally be used in decayers.
   * @param in The particle data objects for the external particles
   * @param intin A pointer to the DecayIntegrator class using this mode.
   * @param onShell Whether or not to perform tests for zero width on-shell particles
   */
  DecayPhaseSpaceMode(tPDVector in, tcDecayIntegratorPtr intin,bool onShell=false) 
    :  _integrator(intin), _maxweight(0.),
       _niter(10), _npoint(10000), _ntry(500),
       _extpart(in),  _partial(-1), _testOnShell(onShell), _ichannel(999) {}
  //@}

  /**
   * Access to the external particles.
   * @param ix The external particle required.
   * @return A pointer to the ParticleData object.
   */
  tcPDPtr externalParticles(int ix) const {return _extpart[ix];}

  /**
   * Number of external particles.
   * @return The number of external particles.
   */
  unsigned int numberofParticles() const {return _extpart.size();}

  /**
   * Number of channels
   * @return The number of channels.
   */
  unsigned int numberChannels() const {return _channels.size();}

  /**
   * Add a new channel. 
   * @param channel A pointer to the new DecayPhaseChannel
   */
  void addChannel(DecayPhaseSpaceChannelPtr channel) {
    channel->init();
    _channels.push_back(channel);
  }

  /**
   * Reset the properties of one of the intermediate particles. Only a specific channel
   * is reset.
   * @param ichan The channel to reset.
   * @param part The ParticleData object of the particle to reset
   * @param mass The mass of the intermediate.
   * @param width The width of gthe intermediate.
   */
  void resetIntermediate(int ichan, tcPDPtr part, Energy mass, Energy width) {
    if(!part) return;
    _channels[ichan]->resetIntermediate(part,mass,width);
  }

  /**
   * Reset the properties of one of the intermediate particles. All the channels 
   * are reset.
   * @param part The ParticleData object of the particle to reset
   * @param mass The mass of the intermediate.
   * @param width The width of gthe intermediate.
   */
  void resetIntermediate(tcPDPtr part, Energy mass, Energy width) {
    for(unsigned int ix=0,N=_channels.size();ix<N;++ix)
      _channels[ix]->resetIntermediate(part,mass,width);
  }

  /**
   * Get the maximum weight for the decay.
   * @return The maximum weight.
   */
  double maxWeight() const {return _maxweight;}

  /**
   * Set the maximum weight for the decay.
   * @param wgt The maximum weight. 
   */
  void setMaxWeight(double wgt) const {_maxweight=wgt;}

  /**
   * Get the weight for a channel. This is the weight for the multi-channel approach.
   * @param ichan The channel.
   * @return The weight for the channel.
   */
  double channelWeight(unsigned int ichan) const {return _channelwgts[ichan];}

  /**
   * Set the weights for the different channels.
   * @param in The weights for the different channels in the multi-channel approach.
   */
  void setWeights(const vector<double> & in) const {_channelwgts=in;}

  /**
   *  Access to the selected channel
   */
  unsigned int selectedChannel() const {return _ichannel;}

  /**
   *  test on/off-shell kinematics
   */
  bool testOnShell() const { return _testOnShell; }

  /**
   *  Access to the epsilon parameter
   */
  Energy epsilonPS() const {return _eps;}

protected:

  /** @name Set-up, Initialization and Access Members */
  //@{
  /**
   * Initialise the phase space.
   * @param init Perform the initialization.
   */
  Energy initializePhaseSpace(bool init, bool onShell=false);

  /**
   * Set the integration parameters
   * @param iter The number of iterations to use for initialization.
   * @param points The number of points to use for each iteration during initialization.
   * @param ntry The number of tries to generate a decay.
   */
  void setIntegrate(int iter,int points,int ntry) {
    _niter=iter;
    _npoint=points;
    _ntry=ntry;
  }

  /**
   * Set the partial width to use for normalization. This is the partial width
   * in the WidthGenerator object.
   * @param in The partial width to use.
   */
  void setPartialWidth(int in) {_partial=in;}
  //@}


  /** @name Phase-Space Generation Members */
  //@{

  /**
   * Access to the matrix element from the decayer.
   * @param ichan The channel, this is to allow the matrix element to be used to
   *              select the intermediates
   * @param inpart The incoming particle.
   * @param opt The option for what to calculate
   * @param outpart The outgoing particles.
   */
  double me2(const int ichan ,const Particle & inpart,
	     const ParticleVector &outpart,DecayIntegrator::MEOption opt) const {
    return _integrator->me2(ichan,inpart,outpart,opt);
  }
  
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
  int selectChannel(const Particle & inpart, ParticleVector & outpart) const {
    // if using flat phase-space don't need to do this
    if(_channelwgts.empty()) return 0;
    vector<double> mewgts(_channels.size(),0.0);
    double total=0.;
    for(unsigned int ix=0,N=_channels.size();ix<N;++ix) {
      mewgts[ix]=me2(ix,inpart,outpart,DecayIntegrator::Calculate);
      total+=mewgts[ix];
    }
    // randomly pick a channel
    total*=UseRandom::rnd();
    int ichan=-1;
    do {
      ++ichan;
      total-=mewgts[ichan];
    }
    while(ichan<int(_channels.size())&&total>0.);
    return ichan;
  }

  /**
   * Return the weight for a given phase-space point.
   * @param cc Whether we are generating the mode specified or the charge 
   *           conjugate mode.
   * @param ichan The channel to generate the weight for.
   * @param in The incoming particle.
   * @param particles The outgoing particles.
   * @param first Whether or not this is the first call for initialisation
   * @return The weight.
   */
  Energy weight(bool cc,int & ichan, const Particle & in,
		ParticleVector & particles,bool first,
		bool onShell=false) const {
    ichan=0;
    Energy phwgt = (_channels.size()==0) ? 
      flatPhaseSpace(cc,in,particles,onShell) : channelPhaseSpace(cc,ichan,in,particles,onShell);
    // generate the matrix element
    return me2(-1,in,particles,
	       first ? DecayIntegrator::Initialize : DecayIntegrator::Calculate)*phwgt;
  }
    
  /**
   * Return the weight and momenta for a flat phase-space decay.
   * @param cc Whether we are generating the mode specified or the charge 
   *           conjugate mode.
   * @param inpart The incoming particle.
   * @param outpart The outgoing particles.
   * @return The weight.
   */
  Energy flatPhaseSpace(bool cc,const Particle & inpart, ParticleVector & outpart,
			bool onShell=false) const;
  
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
			   ParticleVector & particles,
			   bool onShell=false) const;

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
  vector<Energy> externalMasses(Energy inmass,double & wgt, bool onShell) const;
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void init();

  /**
   * Initialize this object to the begining of the run phase.
   */
  virtual void initrun();
  //@}

private:

  /**
   * Private and non-existent assignment operator.
   */
  DecayPhaseSpaceMode & operator=(const DecayPhaseSpaceMode &) = delete;

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
  mutable double _maxweight;

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
  tPDVector _extpart;

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

  /**
   *  Whether to check on-shell or off-shell kinematics
   * in doinit, if on-shell off-shell is tested in initrun
   */
  bool _testOnShell;

  /**
   *  The selected channel
   */
  mutable unsigned int _ichannel;

  /**
   *   Epsilon parameter for phase-space integration
   */
  Energy _eps;

};

  /**
   *  The output operator which is used for debugging.
   */
ostream & operator<<(ostream &, const DecayPhaseSpaceMode &);

}


#endif /* HERWIG_DecayPhaseSpaceMode_H */
