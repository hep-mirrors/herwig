// -*- C++ -*-
#ifndef HERWIG_TensorMesonVectorVectorDecayer_H
#define HERWIG_TensorMesonVectorVectorDecayer_H
//
// This is the declaration of the <!id>TensorMesonVectorVectorDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>TensorMesonVectorVectorDecayer<!!id> class is designed to simulate
//  the decay of a tensor meson to two spin-1 particles. In practice, at least
//  for the light tensor mesons, this is only the decay of a tensor meson to two
//  photons. In principle for bottom and charm tensors this may be the decay to
//  two vector mesons.
//
//  The form of the matrix element is based on the perturbative matrix element
//  for the decay of a graviton to two vector bosons with the neglect of a mass term
//
//  M = T^{ab} ( (e^1_c p^1_a - e^1_a p^1_c)(e^2_c p^2_b - e^2_b p^2_c)
//              +(e^1_c p^1_b - e^1_b p^1_c)(e^2_c p^2_a - e^2_b p^2_a)
//              -0.5 g_{ab}(e^1_c p^1_d- e^1_d p^1_c)(e^2_c p^2_d - e^2_d p^2_c))
//
//  in such a way that it vanishes if the polarizations of the outgoing vectors are
//  replaced with their momenta.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="TensorMesonDecayerBase.html">TensorMesonDecayerBase.h</a>.
// 

#include "TensorMesonDecayerBase.h"
// #include "TensorMesonVectorVectorDecayer.fh"
// #include "TensorMesonVectorVectorDecayer.xh"

namespace Herwig {
using namespace ThePEG; 

class TensorMesonVectorVectorDecayer: public TensorMesonDecayerBase {

public: 

  vector<LorentzTensor> decayTensor(const bool, const int, const int, 
				    const Particle &,
				    const ParticleVector &) const;
  // the hadronic tensor

public:

  inline TensorMesonVectorVectorDecayer();
  inline TensorMesonVectorVectorDecayer(const TensorMesonVectorVectorDecayer &);
  virtual ~TensorMesonVectorVectorDecayer();
  // Standard ctors and dtor.

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

  static ClassDescription<TensorMesonVectorVectorDecayer> initTensorMesonVectorVectorDecayer;
  // Describe a concrete class with persistent data.

  TensorMesonVectorVectorDecayer & operator=(const TensorMesonVectorVectorDecayer &);
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
// base class of Herwig::TensorMesonVectorVectorDecayer.
template <>
struct BaseClassTrait<Herwig::TensorMesonVectorVectorDecayer,1> {
  typedef Herwig::TensorMesonDecayerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::TensorMesonVectorVectorDecayer>
  : public ClassTraitsBase<Herwig::TensorMesonVectorVectorDecayer> {
  static string className() { return "/Herwig++/TensorMesonVectorVectorDecayer"; }
  // Return the class name.
  static string library() { return "libHwTMDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "TensorMesonVectorVectorDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMesonVectorVectorDecayer.tcc"
#endif

#endif /* HERWIG_TensorMesonVectorVectorDecayer_H */
