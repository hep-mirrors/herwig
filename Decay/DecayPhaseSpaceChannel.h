// -*- C++ -*-
#ifndef HERWIG_DecayPhaseSpaceChannel_H
#define HERWIG_DecayPhaseSpaceChannel_H
//
// This is the declaration of the <!id>DecayPhaseSpaceChannel<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is designed to store the information needed for a given
// phase-space channel for use by the multi-channel phase space decayer
// and perform the generation of the phase-space for that channel.
//
// The decay channel is specified as a series of 1->2 decays to either
// external particles or other intermediates. For each intermediate
// the Jacobian to be used can be either a Breit-Wigner(0) or a power-law
// (1).
//
// The class is then capable of generating a phase-space point using this
// channel and computing the weight of a given point for use in a multi-channel
// phase-space integration using the <!id>DecayIntegrator<!!id> class.
//
// The class is designed so that the phase-space channels can either by specified
// using the addIntermediate and setExternal methods directly or via the repository.
// (In practice at the moment all the channels are constructed by the relevant decayers
//  using the former method at the moment.)
//
// CLASSDOC SUBSECTION See also:
//
// <a href="DecayIntegrator.html">DecayIntegrator.h</a>.
// 
// Author: Peter Richardson

#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/Repository/CurrentGenerator.h>

namespace Herwig {
using namespace ThePEG;

class DecayPhaseSpaceChannel: public Interfaced {

friend ostream & operator<<(ostream &, const DecayPhaseSpaceChannel &);

public:
    
  inline DecayPhaseSpaceChannel();
  inline DecayPhaseSpaceChannel(const DecayPhaseSpaceChannel &);
  virtual ~DecayPhaseSpaceChannel();
  // Standard ctors and dtor.
  
public:
    
  vector<Lorentz5Momentum> generateMomenta(const Lorentz5Momentum &);
  // generate the momenta of the external particles
  
  double generateWeight(const vector<Lorentz5Momentum> &);
  // generate the weight for this channel given a phase space configuration
  
  inline PDVector externalParticles() const;
  // return the external particles
  
  inline void addIntermediate(PDPtr,int,double,int,int);
  // add a new intermediate particle
  
  inline void setExternal(PDVector);
  // set the external particles
  
  inline void resetIntermediate(tcPDPtr,Energy,Energy);
  // reset the properties of an intermediate particle
  
  void generateIntermediates(const Particle &, ParticleVector &);
  // generate the final-state particles including the intermediate resonances
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
  
protected:
  
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.
  
protected:
  
  inline virtual void doupdate() throw(UpdateException);
  virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.
  
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.
  
  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.
  
private:
  
  static ClassDescription<DecayPhaseSpaceChannel> initDecayPhaseSpaceChannel;
  // Describe a concrete class with persistent data.
  
  DecayPhaseSpaceChannel & operator=(const DecayPhaseSpaceChannel &);
  // Private and non-existent assignment operator.
  
private:
  
  inline Energy generateMass(int,Energy,Energy);
  // generate the mass of a resonance
  
  inline double massWeight(int,Energy,Energy,Energy);
  // return the weight for a given resonance
  
private:
  
  // properties of intermediate particles
  vector <PDPtr> _intpart;
  // intermediate particles
  vector <int> _jactype;
  // type of jacobian
  vector<Energy> _intmass;
  // mass of particle
  vector<Energy> _intwidth;
  // width of particle
  vector<double> _intpower;
  // power for jacobian if needed
  vector<int> _intdau1; 
  //first daughter
  vector<int> _intdau2;
  // second daughter
  vector< vector<int> > _intext;
  // external particles for the intermediates
  
  // properties of the external particles
  vector <PDPtr> _extpart;
  vector <Energy> _extmass;
  
};

ostream & operator<<(ostream &, const DecayPhaseSpaceChannel &);
// write the phase space channel to a stream
}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of DecayPhaseSpaceChannel.
template <>
struct BaseClassTrait<Herwig::DecayPhaseSpaceChannel,1> {
  typedef Interfaced NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::DecayPhaseSpaceChannel>
  : public ClassTraitsBase<Herwig::DecayPhaseSpaceChannel> {
  static string className() { return "/Herwig++/DecayPhaseSpaceChannel"; }
  // Return the class name.
  static string library() { return "libHwnewDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "DecayPhaseSpaceChannel.icc"

#endif /* HERWIG_DecayPhaseSpaceChannel_H */
