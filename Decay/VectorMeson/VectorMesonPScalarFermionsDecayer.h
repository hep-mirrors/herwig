// -*- C++ -*-
#ifndef HERWIG_VectorMesonPScalarFermionsDecayer_H
#define HERWIG_VectorMesonPScalarFermionsDecayer_H
//
// This is the declaration of the <!id>VectorMesonPScalarFermionsDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// The <!id>VectorMesonPScalarFermionsDecayer<!!id> class is designed to perform the
// decay of a vector meson to a pesudo scalar and a fermion-antifermion pair according
// to a current which is the V->VP vertex combined with the branching of the vector
// into a fermion-antifermion pair.
//
//  It includes the option of a VMD type form-factor of M^2/(m^2_{ff}-M^2+iGamma*M)
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VectorMesonDecayerBase.html">VectorMesonDecayerBase.h</a>.
// 
//  Author: Peter Richardson
//

#include "VectorMesonDecayerBase.h"
// #include "VectorMesonPScalarFermionsDecayer.fh"
// #include "VectorMesonPScalarFermionsDecayer.xh"

namespace Herwig {

using namespace ThePEG;

class VectorMesonPScalarFermionsDecayer: public VectorMesonDecayerBase {

public:

  inline VectorMesonPScalarFermionsDecayer();
  inline VectorMesonPScalarFermionsDecayer(const VectorMesonPScalarFermionsDecayer &);
  virtual ~VectorMesonPScalarFermionsDecayer();
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

  static ClassDescription<VectorMesonPScalarFermionsDecayer> initVectorMesonPScalarFermionsDecayer;
  // Describe a concrete class with persistent data.

  VectorMesonPScalarFermionsDecayer & operator=(const VectorMesonPScalarFermionsDecayer &);
  // Private and non-existent assignment operator.

private:

  
  vector<double> _coupling;
  // coupling for a decay

  vector<int> _incoming,_outgoingP,_outgoingf,_outgoinga;
  // the PDG codes for the incoming and outgoing particles

  vector<double> _maxweight;
  // maximum weight for a decay

  vector<int> _includeVMD;
  vector<int> _VMDid;
  vector<Energy> _VMDmass;
  vector<Energy> _VMDwidth;
};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of Herwig::VectorMesonPScalarFermionsDecayer.
template <>
struct BaseClassTrait<Herwig::VectorMesonPScalarFermionsDecayer,1> {
  typedef Herwig::VectorMesonDecayerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::VectorMesonPScalarFermionsDecayer>
  : public ClassTraitsBase<Herwig::VectorMesonPScalarFermionsDecayer> {
  static string className() { return "/Herwig++/VectorMesonPScalarFermionsDecayer"; }
  // Return the class name.
  static string library() { return "libHwVNDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "VectorMesonPScalarFermionsDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonPScalarFermionsDecayer.tcc"
#endif

#endif /* HERWIG_VectorMesonPScalarFermionsDecayer_H */
