// -*- C++ -*-
#ifndef HERWIG_TensorMesonDecayerBase_H
#define HERWIG_TensorMesonDecayerBase_H
//
// This is the declaration of the <!id>TensorMesonDecayerBase<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>TensorMesonDecayerBase<!!id> is designed to be the
//  base class for the decay of TensorMesons. All decayers implementation
//  such decays should inherit from it and implement the virtual decayTensor
//  method.
//
//  It is baqsed on the DecayIntegrator class to provide multichannel 
//  phase-space integration.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="DecayIntegrator.html">decayIntegrator.h</a>,
// <a href="http:.html">.h</a>.
// 

#include "Herwig++/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/LorentzTensor.h"
// #include "TensorMesonDecayerBase.fh"
// #include "TensorMesonDecayerBase.xh"

namespace Herwig {
using ThePEG::Helicity::LorentzTensor;

class TensorMesonDecayerBase: public DecayIntegrator {

public:

  inline TensorMesonDecayerBase();
  inline TensorMesonDecayerBase(const TensorMesonDecayerBase &);
  virtual ~TensorMesonDecayerBase();
  // Standard ctors and dtor.

public:

  virtual vector<LorentzTensor> decayTensor(const bool, const int, const int, 
					    const Particle &,
					    const ParticleVector &) const;
  // the hadronic tensor

  double me2(bool,const int, const int,
	     const Particle &, const ParticleVector &) const;
  // calculation of the matrix element

public:

  virtual bool accept(const DecayMode &) const;
  // return true if this decayer can perfom the decay specified by the
  // given decay mode.

  virtual ParticleVector decay(const DecayMode &, const Particle &) const;
  // for a given decay mode and a given particle instance, perform the
  // decay and return the decay products.


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

  static AbstractClassDescription<TensorMesonDecayerBase> initTensorMesonDecayerBase;
  // Describe an abstract base class with persistent data.

  TensorMesonDecayerBase & operator=(const TensorMesonDecayerBase &);
  // Private and non-existent assignment operator.

};

  // exception to be thrown if an error
  class TensorDecayerError: public Exception {};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of TensorMesonDecayerBase.
template <>
struct BaseClassTrait<Herwig::TensorMesonDecayerBase,1> {
  typedef Herwig::DecayIntegrator NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::TensorMesonDecayerBase>
  : public ClassTraitsBase<Herwig::TensorMesonDecayerBase> {
  static string className() { return "/Herwig++/TensorMesonDecayerBase"; }
  // Return the class name.
  static string library() { return "libHwTMDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "TensorMesonDecayerBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMesonDecayerBase.tcc"
#endif

#endif /* HERWIG_TensorMesonDecayerBase_H */
