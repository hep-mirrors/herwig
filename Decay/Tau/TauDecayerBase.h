// -*- C++ -*-
#ifndef HERWIG_TauDecayerBase_H
#define HERWIG_TauDecayerBase_H
//
// This is the declaration of the <!id>TauDecayerBase<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class is desgined to be the base class for the HERWIG++ tau decays.
//  It handles the generation of the phase space and calculation of the matrix
//  element. All actual implementations should inherit from this class and
//  implement the virtual hadronCurrent method to return the hadronic
//  current for a given phase space point.
//
//  It also includes utilities to calculate the Breit-Wigners for some
//  intermediate particles which commonly occur in the hadronic currents
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:.html">.h</a>,
// <a href="http:.html">.h</a>.
// 
//  Author: Peter Richardson
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/LorentzSpinor.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;

class TauDecayerBase: public DecayIntegrator {
    
public:
  
  inline TauDecayerBase();
  inline TauDecayerBase(const TauDecayerBase &);
  virtual ~TauDecayerBase();
  // Standard ctors and dtor.
  
public:
  
  virtual bool accept(const DecayMode &) const;
  // return true if this decayer can perfom the decay specified by the
  // given decay mode.
  
  virtual ParticleVector decay(const DecayMode &, const Particle &) const;
  // for a given decay mode and a given particle instance, perform the
  // decay and return the decay products.
  
public:
  
  vector<LorentzPolarizationVector> 
  leptonCurrent(bool vertex, const Particle &, const ParticleVector &) const;
  // the lepton currents for the different lepton helicities
  
  virtual vector<LorentzPolarizationVector> 
  hadronCurrent(bool vertex, const int, const int, 
		const Particle &, const ParticleVector &) const;
  // the hadronic currents    
  
  double me2(bool,const int, const int,
	     const Particle &, const ParticleVector &) const;
  // combine the currents to give the matrix element
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
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
  
private:
  
  static AbstractClassDescription<TauDecayerBase> initTauDecayerBase;
  // Describe an abstract base class with persistent data.
  
  TauDecayerBase & operator=(const TauDecayerBase &);
  // Private and non-existent assignment operator.
  
private:
  
  InvEnergy2 _GF;
  // Fermi coupling constant
  
};
  
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of TauDecayerBase.
  template <>
  struct BaseClassTrait<Herwig::TauDecayerBase,1> {
    typedef Herwig::DecayIntegrator NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::TauDecayerBase>
    : public ClassTraitsBase<Herwig::TauDecayerBase> {
    static string className() { return "/Herwig++/TauDecayerBase"; }
    // Return the class name.
    static string library() { return "libHwTauDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "TauDecayerBase.icc"

#endif /* HERWIG_TauDecayerBase_H */
