// -*- C++ -*-
#ifndef HERWIG_VectorMeson2FermionDecayer_H
#define HERWIG_VectorMeson2FermionDecayer_H
//
// This is the declaration of the <!id>VectorMeson2FermionDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>VectorMeson2FermionDecayer<!!id> class is designed for the decay
//  of a vector meson to a fermion-antifermion pair. This class is mainly used for the
//  decay of the vector mesons to electron and muon pairs.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VectorMesonDecayerBase.html">VectorMesonDecayerBase.h</a>.
// 

#include "VectorMesonDecayerBase.h"

namespace Herwig {
using namespace ThePEG;

class VectorMeson2FermionDecayer: public VectorMesonDecayerBase {

public:

  inline VectorMeson2FermionDecayer();
  inline VectorMeson2FermionDecayer(const VectorMeson2FermionDecayer &);
  virtual ~VectorMeson2FermionDecayer();
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

public:
  
  virtual vector<LorentzPolarizationVector> 
  decayCurrent(const bool, const int, const int, 
	       const Particle &, const ParticleVector &) const;
  // the hadronic currents    
  
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

  static ClassDescription<VectorMeson2FermionDecayer> initVectorMeson2FermionDecayer;
  // Describe a concrete class with persistent data.

  VectorMeson2FermionDecayer & operator=(const VectorMeson2FermionDecayer &);
  // Private and non-existent assignment operator.

private:

  vector<double> _coupling;
  // coupling for a decay

  vector<int> _incoming,_outgoing1,_outgoing2;
  // the PDG codes for the incoming and outgoing particles

  vector<double> _maxweight;
  // maximum weight for a decay

};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of VectorMeson2FermionDecayer.
template <>
struct BaseClassTrait<Herwig::VectorMeson2FermionDecayer,1> {
  typedef Herwig::VectorMesonDecayerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::VectorMeson2FermionDecayer>
  : public ClassTraitsBase<Herwig::VectorMeson2FermionDecayer> {
  static string className() { return "/Herwig++/VectorMeson2FermionDecayer"; }
  // Return the class name.
  static string library() { return "libHwVMDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "VectorMeson2FermionDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMeson2FermionDecayer.tcc"
#endif

#endif /* HERWIG_VectorMeson2FermionDecayer_H */
