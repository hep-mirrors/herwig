// -*- C++ -*-
#ifndef HERWIG_ExampleDecayer_H
#define HERWIG_ExampleDecayer_H
//
// This is the declaration of the <!id>ExampleDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// Dummy class which provides just an example of how to inherits
// from the abstract class  Decayer  and then implement the
// physics of the decay in the two methods:
//           accept
//           decay
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:.html">.h</a>,
// <a href="http:.html">.h</a>.
// 

#include <ThePEG/PDT/Decayer.h>


namespace Herwig {


using namespace ThePEG;

class ExampleDecayer: public ThePEG::Decayer {

public:

  inline ExampleDecayer();
  inline ExampleDecayer(const ExampleDecayer &);
  virtual ~ExampleDecayer();
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
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<ExampleDecayer> initExampleDecayer;
  // Describe a concrete class with persistent data.

  ExampleDecayer & operator=(const ExampleDecayer &);
  //  Private and non-existent assignment operator.

};


}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of ExampleDecayer.
template <>
struct BaseClassTrait<Herwig::ExampleDecayer,1> {
  typedef ThePEG::Decayer NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ExampleDecayer>: public ClassTraitsBase<Herwig::ExampleDecayer> {
  static string className() { return "/Herwig++/ExampleDecayer"; }
  // Return the class name.
  static string library() { return "libHwDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "ExampleDecayer.icc"

#endif /* HERWIG_ExampleDecayer_H */
