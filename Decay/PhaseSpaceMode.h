// -*- C++ -*-
#ifndef Herwig_PhaseSpaceMode_H
#define Herwig_PhaseSpaceMode_H
//
// This is the declaration of the PhaseSpaceMode class.
//

#include "ThePEG/Config/ThePEG.h"
#include "PhaseSpaceMode.fh"
#include "PhaseSpaceChannel.h"
#include "Herwig/PDT/GenericWidthGenerator.h"
#include "Herwig/PDT/GenericMassGenerator.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>PhaseSpaceMode</code> class is designed to store a group
 * of phase-space channels for use by the DecayIntegrator class to 
 * generate the phase-space for a given decay mode.
 *
 * Additional phase-space channels can be added using the addChannel member.
 * 
 *  In practice the modes are usually constructed together with the a number of
 *  <code>PhaseSpaceChannel</code> objects. In classes inheriting from the
 *  DecayIntegrator class.
 *
 * @see DecayIntegrator
 * @see PhaseSpaceChannel
 *
 * @author  Peter Richardson
 * 
 */
class PhaseSpaceMode: public Base {

  friend class PhaseSpaceChannel;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PhaseSpaceMode() : maxWeight_(0.), partial_(-1),
		     testOnShell_(false), eMax_(ZERO)
  {};
  
  /**
   * The default constructor.
   */
  PhaseSpaceMode(tPDPtr in1, tPDVector out,
		 double wMax, tPDPtr in2=tPDPtr(),
		 Energy eMax=ZERO) : incoming_(make_pair(in1,in2)),
				     maxWeight_(wMax),
				     outgoing_(out), partial_(-1),
				     testOnShell_(false), eMax_(eMax)
  {};
  //@}

public:
  
  /**
   * Generate the decay.
   * @param intermediates Whether or not to generate the intermediate particle
   *                      in the decay channel.
   * @param cc Whether we are generating the mode specified or the charge 
   *           conjugate mode.
   * @param inpart The incoming particle.
   * @return The outgoing particles.
   */
  ParticleVector generateDecay(const Particle & inpart,
			       tcDecayIntegrator2Ptr decayer,
			       bool intermediates,bool cc);
  
  /**
   * Add a new channel. 
   * @param channel A pointer to the new PhaseChannel
   */
  void addChannel(PhaseSpaceChannel channel) {
    channel.init(this);
    channels_.push_back(channel);
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
    channels_[ichan].resetIntermediate(part,mass,width);
  }

  /**
   * Reset the properties of one of the intermediate particles. All the channels 
   * are reset.
   * @param part The ParticleData object of the particle to reset
   * @param mass The mass of the intermediate.
   * @param width The width of gthe intermediate.
   */
  void resetIntermediate(tcPDPtr part, Energy mass, Energy width) {
    for(PhaseSpaceChannel & channel : channels_)
      channel.resetIntermediate(part,mass,width);
  }

  /**
   *   The phase-space channels
   */
  const vector<PhaseSpaceChannel> & channels() const {return channels_;}

  /**
   *   Set the weights
   */
  void setWeights(const vector<double> & wgts) {
    assert(wgts.size()==channels_.size());
    for(unsigned int ix=0;ix<channels_.size();++ix)
      channels_[ix].weight(wgts[ix]);
  }
  
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

public :

  /**
   *   Initialize the phase space
   */
  void init() {
    // get mass and width generators
    outgoingCC_.clear();
    for(tPDPtr part : outgoing_) {
      outgoingCC_.push_back(part->CC() ? part->CC() : part);
    }
    massGen_.resize(outgoing_.size());
    widthGen_ = dynamic_ptr_cast<cGenericWidthGeneratorPtr>(incoming_.first->widthGenerator());
    for(unsigned int ix=0;ix<outgoing_.size();++ix) {
      assert(outgoing_[ix]);
      massGen_[ix]= dynamic_ptr_cast<cGenericMassGeneratorPtr>(outgoing_[ix]->massGenerator());
    }
    // get max energy for decays
    if(!incoming_.second)
      eMax_ = testOnShell_ ? incoming_.first->mass() : incoming_.first->massMax();
  }

  void initrun() {
    // update the mass and width generators
    if(incoming_.first->widthGenerator()!=widthGen_)
      widthGen_ = dynamic_ptr_cast<cGenericWidthGeneratorPtr>(incoming_.first->widthGenerator());
    for(unsigned int ix=0;ix<outgoing_.size();++ix) {
      if(massGen_[ix]!=outgoing_[ix]->massGenerator())
	massGen_[ix] = dynamic_ptr_cast<cGenericMassGeneratorPtr>(outgoing_[ix]->massGenerator());
    }
    for(PhaseSpaceChannel & channel : channels_) channel.initrun(this);
    if(widthGen_) const_ptr_cast<GenericWidthGeneratorPtr>(widthGen_)->initrun();
    tcGenericWidthGeneratorPtr wtemp;
    for(unsigned int ix=0;ix<outgoing_.size();++ix) {
      wtemp=
	dynamic_ptr_cast<tcGenericWidthGeneratorPtr>(outgoing_[ix]->widthGenerator());
      if(wtemp) const_ptr_cast<tGenericWidthGeneratorPtr>(wtemp)->initrun();
    }
  }

  /**
   * Get the maximum weight for the decay.
   * @return The maximum weight.
   */
  double maxWeight() const {return maxWeight_;}
  
  /**
   * Set the maximum weight for the decay.
   * @return The maximum weight.
   */
  void maxWeight(double wgt) const {maxWeight_=wgt;}
  
  /**
   * Initialise the phase space.
   * @param init Perform the initialization.
   */
  Energy initializePhaseSpace(bool init, tcDecayIntegrator2Ptr decayer,
			      bool onShell=false);

  /**
   *   The incoming particles
   */
  pair<PDPtr,PDPtr> incoming() const {return incoming_;}
  
  /**
   * Access to the outging particles.
   * @return A pointer to the ParticleData object.
   */
  tPDVector outgoing() const {return outgoing_;}

  /**
   * Number of outgoing particles.
   * @return The number of outgoing particles.
   */
  unsigned int numberOfParticles() const {return outgoing_.size();}

  /**
   * Set the partial width to use for normalization. This is the partial width
   * in the WidthGenerator object.
   * @param in The partial width to use.
   */
  void setPartialWidth(int in) {partial_=in;}
  
public :

  /**
   * A friend operator to allow the mode to be outputted for debugging purposes.
   */
  friend ostream & operator<<(ostream & os, const PhaseSpaceMode & mode) {
    os << "The mode has " << mode.channels_.size() << " channels\n";
    if(mode.incoming_.second==PDPtr())
      os << "This is a mode for the decay of " << mode.incoming_.first->PDGName() << " to ";
    else
      os << "This is a mode for " << mode.incoming_.first->PDGName() << ", "
	 << mode.incoming_.second->PDGName() << " to ";
    for(tPDPtr out : mode.outgoing_) os << out->PDGName() << " ";
    os << "\n";
    for(const PhaseSpaceChannel & channel : mode.channels_) os << channel;
    return os;
  }

private: 

  /**
   * Return the weight for a given phase-space point.
   * @param ichan The channel to generate the weight for.
   * @param in The incoming particle.
   * @param particles The outgoing particles.
   * @param first Whether or not this is the first call for initialisation
   * @return The weight.
   */
  Energy weight(int & ichan, const Particle & in,
		vector<Lorentz5Momentum> & momenta,
		bool onShell=false) const {
    ichan=0;
    // flat phase-space
    if(channels_.empty())
      return flatPhaseSpace(in,momenta,onShell);
    // multi-channel
    else
      return channelPhaseSpace(ichan,in,momenta,onShell);
  }
    
  /**
   * Return the weight and momenta for a flat phase-space decay.
   * @param inpart The incoming particle.
   * @param outpart The outgoing particles.
   * @return The weight.
   */
  Energy flatPhaseSpace(const Particle & inpart,
			vector<Lorentz5Momentum> & momenta,
			bool onShell=false) const;
  
  /**
   * Generate a phase-space point using multichannel phase space.
   * @param in The incoming particle.
   * @param particles The outgoing particles.
   * @return The weight.
   */
  Energy channelPhaseSpace(int & ichan, const Particle & in, 
			   vector<Lorentz5Momentum> & momenta,
			   bool onShell=false) const;

  /**
   * Generate the masses of the external particles.
   * @param inmass The mass of the decaying particle.
   * @param wgt The weight for the masses.
   * @return The masses.
   */
  vector<Energy> externalMasses(Energy inmass,double & wgt, bool onShell) const;

  /**
   * Construct the vertex for spin corrections
   * @param in The incoming particle.
   * @param out The outgoing particles.
   */
  void constructVertex(const Particle & in, const ParticleVector & out,
		       tcDecayIntegrator2Ptr decayer) const;
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PhaseSpaceMode & operator=(const PhaseSpaceMode &);

private:

  /**
   *  The incoming particles
   */
  pair<PDPtr,PDPtr> incoming_;

  /**
   *  The phase-space channels
   */
  vector<PhaseSpaceChannel> channels_;

  /**
   *  The maximum weight
   */
  mutable double maxWeight_;

  /**
   *  The external particles
   */
  tPDVector outgoing_;
  tPDVector outgoingCC_;

  /**
   * Which of the partial widths of the incoming particle to use
   */
  int partial_;

  /**
   * The width generator for the incoming particle.
   */
  cGenericWidthGeneratorPtr widthGen_;

  /**
   *  The mass generators for the outgoing particles.
   */
  vector<cGenericMassGeneratorPtr> massGen_;

  /**
   *  Whether to check on-shell or off-shell kinematics
   * in doinit, if on-shell off-shell is tested in initrun
   */
  bool testOnShell_;

  /**
   *   The maximum energy for the mode
   */
  Energy eMax_;
};

}

#endif /* Herwig_PhaseSpaceMode_H */
