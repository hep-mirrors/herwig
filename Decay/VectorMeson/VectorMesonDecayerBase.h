// -*- C++ -*-
#ifndef HERWIG_VectorMesonDecayerBase_H
#define HERWIG_VectorMesonDecayerBase_H
//
// This is the declaration of the <!id>VectorMesonDecayerBase<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class is designed to be the base class for the HERWIG++ decays of
//  vector mesons. It handles the generation of the phase space and the calculation
//  of the matrix element. All the implementations should inherit from this class
//  and implement the decayCurrent method to return the current for a given
//  phase space point. This current is then contracted with the polarization
//  vectors for the decaying meson.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="DecayIntegrator.html">DecayIntegrator.h</a>,
// <a href="http:.html">.h</a>.
// 
//  Author: Peter Richardson
//
#include "Herwig++/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;

class VectorMesonDecayerBase: public DecayIntegrator {

public:
  
  inline VectorMesonDecayerBase();
  inline VectorMesonDecayerBase(const VectorMesonDecayerBase &);
  virtual ~VectorMesonDecayerBase();
  // Standard ctors and dtor.
  
public:
  
  virtual bool accept(const DecayMode &) const;
  // return true if this decayer can perfom the decay specified by the
  // given decay mode.
  
  virtual ParticleVector decay(const DecayMode &, const Particle &) const;
  // for a given decay mode and a given particle instance, perform the
  // decay and return the decay products.
  
public:
  
  virtual vector<LorentzPolarizationVector> 
  decayCurrent(const bool, const int, const int, 
	       const Particle &, const ParticleVector &) const;
  // the hadronic currents    
  
  double me2(bool,const int, const int,
	     const Particle &, const ParticleVector &) const;
  // calculation of the matrix element
  
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
  
  static AbstractClassDescription<VectorMesonDecayerBase> initVectorMesonDecayerBase;
  // Describe an abstract base class with persistent data.
  
  VectorMesonDecayerBase & operator=(const VectorMesonDecayerBase &);
  // Private and non-existent assignment operator.
  
};
  
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of VectorMesonDecayerBase.
  template <>
  struct BaseClassTrait<Herwig::VectorMesonDecayerBase,1> {
    typedef Herwig::DecayIntegrator NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::VectorMesonDecayerBase>
    : public ClassTraitsBase<Herwig::VectorMesonDecayerBase> {
    static string className() { return "/Herwig++/VectorMesonDecayerBase"; }
    // Return the class name.
    static string library() { return "libHwVMDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "VectorMesonDecayerBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonDecayerBase.tcc"
#endif

#endif /* HERWIG_VectorMesonDecayerBase_H */
