// -*- C++ -*-
#ifndef HERWIG_VectorMeson2MesonDecayer_H
#define HERWIG_VectorMeson2MesonDecayer_H
//
// This is the declaration of the <!id>VectorMeson2MesonDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class is the implementation for the decay of a vector meson to 
//  two scalar mesons, the classic example is rho -> pi pi, via a current
//  which is the difference of the momenta of the mesons. Obviously the
//  order of the momenta doesn't matter as it will only effect the sign of the
//  matrix element
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:.html">VectorMesonDecayerBase.h</a>.
// 
//  Author: Peter Richardson
//

#include "VectorMesonDecayerBase.h"

namespace Herwig {

using namespace ThePEG;

class VectorMeson2MesonDecayer: public VectorMesonDecayerBase {
  
public:
  
  inline VectorMeson2MesonDecayer();
  inline VectorMeson2MesonDecayer(const VectorMeson2MesonDecayer &);
  virtual ~VectorMeson2MesonDecayer();
  // Standard ctors and dtor.
  
public:
  
  virtual bool accept(const DecayMode &) const;
  // return true if this decayer can perfom the decay specified by the
  // given decay mode.
  
  virtual ParticleVector decay(const DecayMode &, const Particle &) const;
  // for a given decay mode and a given particle instance, perform the
  // decay and return the decay products.
  
  virtual vector<LorentzPolarizationVector> 
  decayCurrent(const bool, const int, const int, 
	       const Particle &, const ParticleVector &) const;
  // the hadronic currents    
  
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
  
  static ClassDescription<VectorMeson2MesonDecayer> initVectorMeson2MesonDecayer;
  // Describe a concrete class with persistent data.
  
  VectorMeson2MesonDecayer & operator=(const VectorMeson2MesonDecayer &);
  // Private and non-existent assignment operator.
  
private:
  
  vector<int> _incoming,_outgoing1,_outgoing2;
  // the PDG codes for the incoming and outgoing particles

  vector<double> _maxweight;
  // the maximum weight for the integration

  vector<double> _coupling;
  // the coupling for the decay

};
  
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of VectorMeson2MesonDecayer.
  template <>
  struct BaseClassTrait<Herwig::VectorMeson2MesonDecayer,1> {
    typedef Herwig::VectorMesonDecayerBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::VectorMeson2MesonDecayer>
    : public ClassTraitsBase<Herwig::VectorMeson2MesonDecayer> {
    static string className() { return "/Herwig++/VectorMeson2MesonDecayer"; }
    // Return the class name.
    static string library() { return "libHwVMDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "VectorMeson2MesonDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMeson2MesonDecayer.tcc"
#endif

#endif /* HERWIG_VectorMeson2MesonDecayer_H */
