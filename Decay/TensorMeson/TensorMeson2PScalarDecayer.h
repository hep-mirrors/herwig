// -*- C++ -*-
#ifndef HERWIG_TensorMeson2PScalarDecayer_H
#define HERWIG_TensorMeson2PScalarDecayer_H
//
// This is the declaration of the <!id>TensorMeson2PScalarDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>TensorMeson2PScalarDecayer<!!id> class is designed for thew decay
//  of a tensor meson to two pseudoscalars via a decay matrix element which 
//  takes the form
// 
//   e^{ab} p^1_a p^2_b
//
//  It can also be used for the decay of a tensor to two scalar mesons although
//  this rarely happens in practice.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="TensorMesonDecayrerBase.html">TensorMesonDecayerBase.h</a>.
// 

#include "TensorMesonDecayerBase.h"
// #include "TensorMeson2PScalarDecayer.fh"
// #include "TensorMeson2PScalarDecayer.xh"

namespace Herwig {
using namespace ThePEG;

class TensorMeson2PScalarDecayer: public TensorMesonDecayerBase {

public:

  inline TensorMeson2PScalarDecayer();
  inline TensorMeson2PScalarDecayer(const TensorMeson2PScalarDecayer &);
  virtual ~TensorMeson2PScalarDecayer();
  // Standard ctors and dtor.

public: 

  vector<LorentzTensor> decayTensor(const bool, const int, const int, 
				    const Particle &,
				    const ParticleVector &) const;
  // the hadronic tensor

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

  static ClassDescription<TensorMeson2PScalarDecayer> initTensorMeson2PScalarDecayer;
  // Describe a concrete class with persistent data.

  TensorMeson2PScalarDecayer & operator=(const TensorMeson2PScalarDecayer &);
  // Private and non-existent assignment operator.

private:

  vector<int> _incoming;
  // the PDG codes for the incoming particles

  vector<int> _outgoing1,_outgoing2;
  // the PDG codes for the outgoing particles

  vector<double> _coupling;
  // the coupling for the decay

  vector<double> _maxweight;
  // the maximum weight for the decay

};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of TensorMeson2PScalarDecayer.
template <>
struct BaseClassTrait<Herwig::TensorMeson2PScalarDecayer,1> {
  typedef Herwig::TensorMesonDecayerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::TensorMeson2PScalarDecayer>
  : public ClassTraitsBase<Herwig::TensorMeson2PScalarDecayer> {
  static string className() { return "/Herwig++/TensorMeson2PScalarDecayer"; }
  // Return the class name.
  static string library() { return "libHwTMDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "TensorMeson2PScalarDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMeson2PScalarDecayer.tcc"
#endif

#endif /* HERWIG_TensorMeson2PScalarDecayer_H */
