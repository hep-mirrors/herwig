// -*- C++ -*-
#ifndef HERWIG_PScalarVectorVectorDecayer_H
#define HERWIG_PScalarVectorVectorDecayer_H
//
// This is the declaration of the <!id>PScalarVectorVectorDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// The <!id>PScalarVectorVectorDecayer<!!id> class is designed to perform the decay 
// of a pseudoscalar meson to two spin-1 particles. The majority of these decays
// are of a light pseudoscalar meson to gamma-gamma where including the matrix-element
// is unnessecary. However there are a small number of decays, eg eta' -> omega gamma
// where the used of this decayer is required to get the correct correlations.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="DecayIntegrator.html">DecayIntegrator.h</a>.
// 

#include "Herwig++/Decay/DecayIntegrator.h"
// #include "PScalarVectorVectorDecayer.fh"
// #include "PScalarVectorVectorDecayer.xh"

namespace Herwig {

class PScalarVectorVectorDecayer: public DecayIntegrator {

public:

  inline PScalarVectorVectorDecayer();
  inline PScalarVectorVectorDecayer(const PScalarVectorVectorDecayer &);
  virtual ~PScalarVectorVectorDecayer();
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

  static ClassDescription<PScalarVectorVectorDecayer> initPScalarVectorVectorDecayer;
  // Describe a concrete class with persistent data.

  PScalarVectorVectorDecayer & operator=(const PScalarVectorVectorDecayer &);
  // Private and non-existent assignment operator.

private:

  vector<int> _incoming;
  // the PDG code for the incoming particle
  vector<int> _outgoing1;
  // the PDG code for the first outgoing particle
  vector<int> _outgoing2;
  // the PDG code for the second incoming particle
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
// base class of PScalarVectorVectorDecayer.
template <>
struct BaseClassTrait<Herwig::PScalarVectorVectorDecayer,1> {
  typedef Herwig::DecayIntegrator NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::PScalarVectorVectorDecayer>
  : public ClassTraitsBase<Herwig::PScalarVectorVectorDecayer> {
  static string className() { return "/Herwig++/PScalarVectorVectorDecayer"; }
  // Return the class name.
  static string library() { return "libHwSMDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "PScalarVectorVectorDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalarVectorVectorDecayer.tcc"
#endif

#endif /* HERWIG_PScalarVectorVectorDecayer_H */
