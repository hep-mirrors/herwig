// -*- C++ -*-
#ifndef HERWIG_TensorMesonVectorPScalarDecayer_H
#define HERWIG_TensorMesonVectorPScalarDecayer_H
//
// This is the declaration of the <!id>TensorMesonVectorPScalarDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>TensorMesonVectorPScalarDecayer<!!id> class handles the decay of
//  a tensor meson to a vector and a pseudoscalar. Examples of this decay
//  are a+2 -> rho pi or a_2 pi gamma.
//
//  The class inherits from the <id>TensorMesonDecayerBase<!id> class and implements
//  the virtual decayTensor method to calculate the tensor for the decay.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="TensorMesonDecayerBase.html">TensorMesonDecayer.h</a>.
// 

#include "TensorMesonDecayerBase.h"
// #include "TensorMesonVectorPScalarDecayer.fh"
// #include "TensorMesonVectorPScalarDecayer.xh"

namespace Herwig {
using namespace ThePEG;

class TensorMesonVectorPScalarDecayer: public TensorMesonDecayerBase {

public:

  inline TensorMesonVectorPScalarDecayer();
  inline TensorMesonVectorPScalarDecayer(const TensorMesonVectorPScalarDecayer &);
  virtual ~TensorMesonVectorPScalarDecayer();
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

  static ClassDescription<TensorMesonVectorPScalarDecayer> initTensorMesonVectorPScalarDecayer;
  // Describe a concrete class with persistent data.

  TensorMesonVectorPScalarDecayer & operator=(const TensorMesonVectorPScalarDecayer &);
  // Private and non-existent assignment operator.

private:

  vector<int> _incoming;
  // PDG codes for the incoming particles

  vector<int> _outgoingV,_outgoingP;
  // PDG codes for the outgoing particles

  vector<double> _coupling;
  // coupling for the decay

  vector<double> _maxweight;
  // max weight ofr the decay
};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of Herwig::TensorMesonVectorPScalarDecayer.
template <>
struct BaseClassTrait<Herwig::TensorMesonVectorPScalarDecayer,1> {
  typedef Herwig::TensorMesonDecayerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::TensorMesonVectorPScalarDecayer>
  : public ClassTraitsBase<Herwig::TensorMesonVectorPScalarDecayer> {
  static string className() { return "/Herwig++/TensorMesonVectorPScalarDecayer"; }
  // Return the class name.
  static string library() { return "libHwTMDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "TensorMesonVectorPScalarDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMesonVectorPScalarDecayer.tcc"
#endif

#endif /* HERWIG_TensorMesonVectorPScalarDecayer_H */
