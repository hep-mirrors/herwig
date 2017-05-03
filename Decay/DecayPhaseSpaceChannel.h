// -*- C++ -*-
//
// DecayPhaseSpaceChannel.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DecayPhaseSpaceChannel_H
#define HERWIG_DecayPhaseSpaceChannel_H
//
// This is the declaration of the DecayPhaseSpaceChannel class.
//
#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/EventRecord/Particle.h>
#include "DecayPhaseSpaceChannel.fh"
#include "DecayIntegrator.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Repository/UseRandom.h"
#include "DecayPhaseSpaceMode.fh"

namespace Herwig {
using namespace ThePEG;

  /** \ingroup Decay
   *
   * This class is designed to store the information needed for a given
   * phase-space channel for use by the multi-channel phase space decayer
   * and perform the generation of the phase space for that channel.
   *
   * The decay channel is specified as a series of \f$1\to2\f$ decays to either
   * external particles or other intermediates. For each intermediate
   * the Jacobian to be used can be either a Breit-Wigner(0) or a power-law
   * (1).
   *
   * The class is then capable of generating a phase-space point using this
   * channel and computing the weight of a given point for use in a multi-channel
   * phase space integration using the <code>DecayPhaseSpaceMode</code> class.
   *
   * The class is designed so that the phase-space channels can either by specified
   * using the addIntermediate method directly or via the repository.
   * (In practice at the moment all the channels are constructed by the relevant decayers
   *  using the former method at the moment.)
   *
   * @see DecayPhaseSpaceMode
   * @see DecayIntegrator
   *
   * @author  Peter Richardson
   */

class DecayPhaseSpaceChannel: public Interfaced {

  /**
   *  A friend output operator to allow the channel to be outputted for
   * debugging purposes
   */
  friend ostream & operator<<(ostream &, const DecayPhaseSpaceChannel &);

  /**
   * DecayPhaseSpaceMode is a friend to avoid making many of the phase space
   * generation members public.
   */
  friend class DecayPhaseSpaceMode;

public:
    
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  DecayPhaseSpaceChannel() {}

  /**
   * Constructor with a pointer to a <code>DecayPhaseSpaceMode</code>. This 
   * is the constructor which should normally be used in decayers.
   */
  DecayPhaseSpaceChannel(tcDecayPhaseSpaceModePtr);
  //@}
  
public:
  
  /** @name Set-up Members */
  //@{
  /**
   * Add a new intermediate particle
   * @param part A pointer to the particle data object for the intermediate.
   * @param jac  The jacobian to be used for the generation of the particle's mass
   * 0 is a Breit-Wigner and 1 is a power-law
   * @param power The power to beb used for the mass generation if a power law
   * mass distribution is chosen.
   * @param dau1 The first daughter. If this is postive it is the \f$dau1\f$ th
   * outgoing particle (0 is the incoming particle), if it is negative it is the 
   * \f$|dau1|\f$ intermediate. The intermediates are specified in the order they
   * are added with 0 being the incoming particle.
   * @param dau2 The first daughter. If this is postive it is the \f$dau2\f$ th
   * outgoing particle (0 is the incoming particle), if it is negative it is the 
   * \f$|dau2|\f$ intermediate. The intermediates are specified in the order they
   * are added with 0 being the incoming particle.
   */
  void addIntermediate(PDPtr part,int jac,double power,int dau1,int dau2) {
    _intpart.push_back(part);
    _jactype.push_back(jac);
    _intpower.push_back(power);
    _intdau1.push_back(dau1); 
    _intdau2.push_back(dau2);
  }
  
  /**
   * Reset the properties of an intermediate particle. This member is used
   * when a Decayer is used a different value for either the mass or width
   * of a resonace to that in the ParticleData object. This improves the 
   * efficiency of the integration.
   * @param part A pointer to the particle data object for the intermediate.
   * @param mass The new mass of the intermediate
   * @param width The new width of the intermediate.
   */
  void resetIntermediate(tcPDPtr part,Energy mass,Energy width) {
    if(!part) return;
    int idin=part->id();
    for(unsigned int ix=0;ix<_intpart.size();++ix) {
      if(_intpart[ix] && _intpart[ix]->id()==idin) {
	_intmass[ix] =mass;_intwidth[ix]=width;
	_intmass2[ix]=mass*mass;_intmwidth[ix]=mass*width;
      }
    }
  }

  /*
   * Reset the one of the daughters
   * @param oldp The id of the particle being reset
   * @param newp The id of the particle replacing it
   */
  void resetDaughter(int oldp, int newp) {
    for(unsigned int ix=0;ix<_intdau1.size();++ix) {
      if(_intdau1[ix]==oldp) _intdau1[ix]=newp;
    }
    for(unsigned int ix=0;ix<_intdau2.size();++ix) {
      if(_intdau2[ix]==oldp) _intdau2[ix]=newp;
    }
  }
  //@}

protected:

  /** @name Phase-Space Generation Members */
  //@{
  /**
   * Generate the momenta of the external particles. This member uses the masses
   * of the external particles generated by the DecayPhaseMode class and the
   * intermediates for the channel to generate the momenta of the decay products.
   * @param pin The momenta of the decay products. This is outputed by the member.
   * @param massext The masses of the particles. This is to allow inclusion of
   * off-shell effects for the external particles.
   */
  vector<Lorentz5Momentum> generateMomenta(const Lorentz5Momentum & pin,
					   const vector<Energy> & massext);
  
  /**
   * Generate the weight for this channel given a phase space configuration.
   * This member generates the weight for a given phase space configuration
   * and is used by the DecayPhaseSpaceMode class to compute the denominator
   * of the weight in the multi-channel integration.
   * @param output The momenta of the outgoing particles.
   */
  double generateWeight(const vector<Lorentz5Momentum> & output);

  /**
   * Generate the final-state particles including the intermediate resonances.
   * This method takes the outgoing particles and adds the intermediate particles
   * specified by this phase-space channel. This is to allow a given set of 
   * intermediates to be specified even when there is interference between different
   * intermediate states.
   * @param cc Whether the particles are the mode specified or its charge conjugate.
   * @param in The incoming particles.
   * @param out The outgoing particles.
   * 
   */
  void generateIntermediates(bool cc,const Particle & in, ParticleVector & out);

  /**
   * Calculate the momenta for a two body decay
   * The return value indicates success or failure.
   * @param p The momentum of the decaying particle
   * @param m1 The mass of the first decay product
   * @param m2 The mass of the second decay product
   * @param p1 The momentum of the first decay product
   * @param p2 The momentum of the second decay product
   */
  void twoBodyDecay(const Lorentz5Momentum & p, 
		    const Energy m1, const Energy m2,
		    Lorentz5Momentum & p1, Lorentz5Momentum & p2);
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
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<DecayPhaseSpaceChannel> initDecayPhaseSpaceChannel;
  
  /**
   * Private and non-existent assignment operator.
   */
  DecayPhaseSpaceChannel & operator=(const DecayPhaseSpaceChannel &);
  
private:
  
  /** @name Mass Generation Members */
  //@{
  /**
   * Generate the mass of a resonance.
   * @param ires The resonance to be generated.
   * @param lower The lower limit on the particle's mass.
   * @param upper The upper limit on the particle's mass. 
   */
  Energy generateMass(int ires,Energy lower,Energy upper);
  
  /**
   * Return the weight for a given resonance.
   * @param ires The resonance to be generated.
   * @param moff The mass of the resonance.
   * @param lower The lower limit on the particle's mass.
   * @param upper The upper limit on the particle's mass. 
   */
  InvEnergy2 massWeight(int ires, Energy moff,Energy lower,Energy upper);
  //@}  

private:
  
  /**
   * pointer to the mode
   */
  tcDecayPhaseSpaceModePtr _mode;

  /**
   * Pointers to the particle data objects of the intermediate particles
   */
  vector <PDPtr> _intpart;

  /**
   * The type of jacobian to be used for the intermediates.
   */
  vector <int> _jactype;

  /**
   * The mass of the intermediates.
   */
  vector<Energy> _intmass;

  /**
   * The width of the intermediates.
   */
  vector<Energy> _intwidth;

  /**
   * The mass squared of the intermediate particles.
   */
  vector<Energy2> _intmass2;

  /**
   * The mass times the width for the intermediate particles.
   */
  vector<Energy2>_intmwidth;

  /**
   * The power for the intermediate resonance.
   */
  vector<double> _intpower;

  /**
   * The first daughter of the intermediate resonance.
   */
  vector<int> _intdau1; 

  /**
   *  The second daughter of the intermediate resonance.
   */
  vector<int> _intdau2;

  /**
   * The external particles that an intermediate particle final decays to.
   */
  vector< vector<int> > _intext;
  
  /**
   * Helper function for the weight calculation.
   * @param ires The resonance to be generated.
   * @param limit The limit on the particle's mass. 
   */
  double atanhelper_(int ires, Energy limit);
};


/**
 * write the phase space channel to a stream
 */
ostream & operator<<(ostream &, const DecayPhaseSpaceChannel &);

/**
 * exception for this class and those inheriting from it
 */
class DecayPhaseSpaceError: public Exception {};
}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of DecayPhaseSpaceChannel.
 */
template <>

struct BaseClassTrait<Herwig::DecayPhaseSpaceChannel,1> {
  /** Typedef of the base class of DecayPhaseSpaceChannel */
  typedef Interfaced NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>

struct ClassTraits<Herwig::DecayPhaseSpaceChannel>
  : public ClassTraitsBase<Herwig::DecayPhaseSpaceChannel> {
    /**  Return the class name.*/
  static string className() { return "Herwig::DecayPhaseSpaceChannel"; }
};

/** @endcond */

}

#endif /* HERWIG_DecayPhaseSpaceChannel_H */
