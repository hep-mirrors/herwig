// -*- C++ -*-
#ifndef HERWIG_DecayIntegrator_H
#define HERWIG_DecayIntegrator_H
//
// This is the declaration of the <!id>DecayIntegrator<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class is designed to be the base class for Herwig++ decays including
//  the implementation of a multichannel decayer or n-body phase space decays.
//
//  The <!id>DecayIntegrator<!!id> class inherits from ThePEG's Decayer class
//  and makes use of the <!id>DecayPhaseChannel<!!id> class to specify the 
//  phase-space channel. 
//
//  Additional phase-space channels can be added using the addChannel member
//  and a group of channels can be specified for a given decay mode, 
//  if a decayer handles more than one mode, using the setMode member.
//  In practice the channels are usually specified in the doinit() member of
//  the Decayer. For the majority of the decays currently implemented the 
//  phase-space integration has been optimised and the maximum weight set.
//  If the parameters of the decay model are changed the initializePhaseSpace
//  member can be called to optimise the integration and calculate the maximum weight.
//  In general this is called in the initrun() member of the decayer for the different
//  decay modes. The phase space initialisation only occurs if the Initialize interface
//  is set to one.
//
//  In classes inheriting for this the me2() member which gives the matrix element
//  squared must be implemented. This should be combined with the setting of the
//  phase space channels, and the setting of which channels to use and their
//  initial weights in the doinit() member. The different decay modes should then
//  be initialized in the initrun() member if needed. The generate member can then
//  be called from the decay() member to generate a phase-space configuration for a 
//  decay.
//   
// CLASSDOC SUBSECTION See also:
//
// <a href="DecayPhaseSpaceChannel.html">.h</a>,
// 
// Author: Peter Richardson

#include <ThePEG/PDT/Decayer.h>
#include "DecayPhaseSpaceChannel.h"
#include <ThePEG/PDT/EnumParticles.h>
#include <Herwig++/Helicity/Correlations/DecayVertex.h>
#include "ThePEG/Utilities/Timer.h"
#include <ThePEG/Helicity/SpinInfo.h>

namespace Herwig {
using namespace ThePEG;
using Herwig::Helicity::DecayMatrixElement;
class DecayIntegrator: public Decayer {
    
friend ostream & operator<<(ostream &, const DecayIntegrator &);

public:
  
  inline DecayIntegrator();
  inline DecayIntegrator(const DecayIntegrator &);
  virtual ~DecayIntegrator();
  // Standard ctors and dtor.
      
public:
  
  virtual bool accept(const DecayMode &) const;
  // return true if this decayer can perfom the decay specified by the
  // given decay mode.
  
  virtual ParticleVector decay(const DecayMode &, const Particle &) const;
  // for a given decay mode and a given particle instance, perform the
  // decay and return the decay products.
  
  inline void setMode(const unsigned int, double, const vector<bool>, 
		      const vector<double>) const;
  // set the parameters for a given mode
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  virtual double me2(bool, const int, const int,
		     const Particle &, const ParticleVector &) const;
  // return the matrix element squared for a given mode and phase-space channel

  inline void constructVertex(const Particle &, const ParticleVector &) const;
  // construct the vertex for spin corrections
    
  double flatPhaseSpace(const Particle &, ParticleVector &) const;
  // return the weight and momenta for a flat phase-space decay
  
  void initializePhaseSpace(const unsigned int imode, const PDVector &);
  // initialise the phase space
  
  inline double getMaxWeight(const unsigned int) const;
  // get the maximum weight for the decay
  
  inline double weight(unsigned int imode, 
		       int &, const Particle &, ParticleVector &) const;
  // return the weight for a given point
  
  void generate(bool,const unsigned int,const Particle &, ParticleVector &) const;
  // generate the decay
  
  double channelPhaseSpace(unsigned int,
			   int &, const Particle &, ParticleVector &) const;
  // generate a phase-space point using multichannel phase space
  
  inline void addChannel(Ptr<Herwig::DecayPhaseSpaceChannel>::pointer);
  // add a new channel
  
  inline unsigned int numberChannels() const;
  // return the number of channels
  
  inline vector<double> getWeights(const unsigned int) const ;
  // get the weights for the different channels
  
  inline void resetIntermediate(int ichan, tcPDPtr part, Energy mass, Energy width);
  // reset the properties of one of the intermediate particles
  
protected:
  
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.
  
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.
  
  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.
  
protected:
  
  inline const DecayMatrixElement & ME() const;
  inline void ME(const DecayMatrixElement &) const;
  // set and get the matrix element
  
private:
  
  inline void setMaxWeight(const unsigned int, double) const;
  // set the maximum weight for the integration channel
  
  inline void setOn(const unsigned int, const vector<bool>) const;
  // set which channels are on
  
  inline void setWeights(const unsigned int,const vector<double>) const;
  // set the weights for the different channels
  
  inline int selectChannel(const unsigned int,
			   const Particle &, ParticleVector &) const;
  // select which channel to use to output the particles

private:
  
  static AbstractClassDescription<DecayIntegrator> initDecayIntegrator;
  // Describe an abstract base class with persistent data.
  
  DecayIntegrator & operator=(const DecayIntegrator &);
  // Private and non-existent assignment operator.
  
private:
  
  vector<Ptr<Herwig::DecayPhaseSpaceChannel>::pointer> _channels;
  // pointers to the phase-space channels
  
  mutable vector< vector<double> > _channelwgts;
  // weights for the channels for the different modes
  mutable vector< vector<bool> > _channelon;
  // which channels are used for the different modes    
  
  mutable vector<double> _MaxWeight;
  // the maximum weight for the integration
  bool _Initialize;
  // perform initialisation
  int _niter,_npoint;
  // parameters for the initialisation
  
  int _ntry;
  // number of attempts to generate a decay
  
  mutable DecayMatrixElement _matrixelement;
  
};
  
ostream & operator<<(ostream &, const DecayIntegrator &);

}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of DecayIntegrator.
  template <>
  struct BaseClassTrait<Herwig::DecayIntegrator,1> {
    typedef Decayer NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::DecayIntegrator>
    : public ClassTraitsBase<Herwig::DecayIntegrator> {
    static string className() { return "/Herwig+/DecayIntegrator"; }
    // Return the class name.
    static string library() { return "libHwnewDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };

}

#include "DecayIntegrator.icc"

#endif /* HERWIG_DecayIntegrator_H */
