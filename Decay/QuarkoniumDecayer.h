// -*- C++ -*-
#ifndef HERWIG_QuarkoniumDecayer_H
#define HERWIG_QuarkoniumDecayer_H
//
// This is the declaration of the <!id>QuarkoniumDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// <!id>QuarkoniumDecayer<!!id> is a class that defines all the general routines 
// used in HERWIG++ to imitate the HERWIG 6.4 decays. The goal is to have an exact
// copy of HERWIG 6.4 decay routines. This will allow for easy 'callibration'
// of the new C++ code with the old Fortran code.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="Decayer.html">Decayer.h</a>.
// 

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/PDT/Decayer.h>
#include <ThePEG/Handlers/HandlerBase.h>
#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Repository/Strategy.fh>
#include "DecayConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include <fstream>

namespace Herwig {

using namespace ThePEG;

class QuarkoniumDecayer: public Decayer {

public:

  inline QuarkoniumDecayer();
  inline QuarkoniumDecayer(const QuarkoniumDecayer &);
  virtual ~QuarkoniumDecayer();
  // Standard ctors and dtor

public:

  virtual bool accept(const DecayMode &) const;
  // return true if this decayer can perfom the decay specified by the
  // given decay mode.

  virtual ParticleVector decay(const DecayMode &, const Particle &) const;
  // for a given decay mode and a given particle instance, perform the
  // decay and return the decay products.

  static void Init();
  // Standard Init function used to initialize the interface.

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  // Standard Interfaced virtual functions

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard Persisten stream methods

protected:
   inline virtual IBPtr clone() const;
   inline virtual IBPtr fullclone() const;
   // Standard clone methods

   inline virtual void rebind(const TranslationMap &trans) 
                              throw(RebindException);
   // Change all pointers of Interfaced objects to corresponding clones

   inline virtual IVector getReferences();
   // Return pointers to all interfaced objects referred to by this class.

private:

  int MECode;
  Ptr<GlobalParameters>::pointer global;

  static ClassDescription<QuarkoniumDecayer> initQuarkoniumDecayer;
  static long lastAddedNumber;

  const QuarkoniumDecayer & operator=(const QuarkoniumDecayer &);
  //  Private and non-existent assignment operator.
};

}

namespace ThePEG {

template <>
struct BaseClassTrait<Herwig::QuarkoniumDecayer,1> {
  typedef HandlerBase NthBase;
};

template <>
struct ClassTraits<Herwig::QuarkoniumDecayer>: public ClassTraitsBase<Herwig::QuarkoniumDecayer> {
  static string className() { return "/Herwig++/QuarkoniumDecayer"; }
  static string library() { return "libHwDecay.so"; }
};

}

#include "QuarkoniumDecayer.icc"

#endif /* HERWIG_QuarkoniumDecayer_H */
