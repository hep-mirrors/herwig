// -*- C++ -*-
#ifndef HERWIG_VectorMesonVectorPScalarDecayer_H
#define HERWIG_VectorMesonVectorPScalarDecayer_H
//
// This is the declaration of the <!id>VectorMesonVectorPScalarDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class is designed for the decay of a vector meson to another spin-1 
//  particle, either another vector meson or a photon, and a pseduoscalar meson.
//  The current for the decay is 
//
//   eps^{a,b,c,d} p0_a eps0_b p1_c eps1_d.
//
//  Examples of such decays are rho -> pi gamma etc,
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VectorMesonDecayerBase.html">VectorMesonDecayerBase.h</a>.
// 
//  Author: Peter Richardson
//

#include "VectorMesonDecayerBase.h"
// #include "VectorMesonVectorPScalarDecayer.fh"
// #include "VectorMesonVectorPScalarDecayer.xh"

namespace Herwig {

class VectorMesonVectorPScalarDecayer: public VectorMesonDecayerBase {

public:

  inline VectorMesonVectorPScalarDecayer();
  inline VectorMesonVectorPScalarDecayer(const VectorMesonVectorPScalarDecayer &);
  virtual ~VectorMesonVectorPScalarDecayer();
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
  // the hadronic current

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

  static ClassDescription<VectorMesonVectorPScalarDecayer> initVectorMesonVectorPScalarDecayer;
  // Describe a concrete class with persistent data.

  VectorMesonVectorPScalarDecayer & operator=(const VectorMesonVectorPScalarDecayer &);
  // Private and non-existent assignment operator.

private:

  vector<double> _coupling;
  // coupling for a decay

  vector<int> _incoming,_outgoingV,_outgoingP;
  // the PDG codes for the incoming and outgoing particles

  vector<double> _maxweight;
  // maximum weight for a decay

};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of VectorMesonVectorPScalarDecayer.
template <>
struct BaseClassTrait<Herwig::VectorMesonVectorPScalarDecayer,1> {
  typedef Herwig::VectorMesonDecayerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::VectorMesonVectorPScalarDecayer>
  : public ClassTraitsBase<Herwig::VectorMesonVectorPScalarDecayer> {
  static string className() { return "/Herwig++/VectorMesonVectorPScalarDecayer"; }
  // Return the class name.
  static string library() { return "libHwVMDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "VectorMesonVectorPScalarDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonVectorPScalarDecayer.tcc"
#endif

#endif /* HERWIG_VectorMesonVectorPScalarDecayer_H */
