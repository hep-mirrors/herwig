// -*- C++ -*-
#ifndef THEPEG_PScalarVectorFermionsDecayer_H
#define THEPEG_PScalarVectorFermionsDecayer_H
//
// This is the declaration of the <!id>PScalarVectorFermionsDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// The <!id>PScalarVectorFermionsDecayer<!!id> class is designed for the decay of a 
// pseudoscalar meson to a spin-1 particle and a fermion-antifermion pair. In practice
// these decays are of the form gamma l+l- and the propagator of the off-shell boson
// is taken to be 1/m^2. There is also the option of including a vector meson dominance
// form-factor.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="DecayIntegrator.html">DecayIntegrato.h</a>.
// 
//  Author: Peter Richardson
//

#include "Herwig++/Decay/DecayIntegrator.h"
// #include "PScalarVectorFermionsDecayer.fh"
// #include "PScalarVectorFermionsDecayer.xh"

namespace Herwig {
using namespace ThePEG;

class PScalarVectorFermionsDecayer: public DecayIntegrator {

public:

  inline PScalarVectorFermionsDecayer();
  inline PScalarVectorFermionsDecayer(const PScalarVectorFermionsDecayer &);
  virtual ~PScalarVectorFermionsDecayer();
  // Standard ctors and dtor.

public:

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

  static ClassDescription<PScalarVectorFermionsDecayer> initPScalarVectorFermionsDecayer;
  // Describe a concrete class with persistent data.

  PScalarVectorFermionsDecayer & operator=(const PScalarVectorFermionsDecayer &);
  // Private and non-existent assignment operator.

private:

  vector<double> _coupling;
  // coupling for a decay

  vector<int> _incoming,_outgoingV,_outgoingf,_outgoinga;
  // the PDG codes for the incoming and outgoing particles

  vector<double> _maxweight;
  // maximum weight for a decay

  vector<int> _includeVMD;
  vector<int> _VMDid;
  vector<Energy> _VMDmass;
  vector<Energy> _VMDwidth;
  // parameters for the VMD factor
};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of PScalarVectorFermionsDecayer.
template <>
struct BaseClassTrait<Herwig::PScalarVectorFermionsDecayer,1> {
  typedef Herwig::DecayIntegrator NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::PScalarVectorFermionsDecayer>
  : public ClassTraitsBase<Herwig::PScalarVectorFermionsDecayer> {
  static string className() { return "/Herwig++/PScalarVectorFermionsDecayer"; }
  // Return the class name.
  static string library() { return "libHwSMDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "PScalarVectorFermionsDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalarVectorFermionsDecayer.tcc"
#endif

#endif /* THEPEG_PScalarVectorFermionsDecayer_H */
