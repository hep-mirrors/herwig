// -*- C++ -*-
#ifndef HERWIG_HeavyDecayer_H
#define HERWIG_HeavyDecayer_H
//
// This is the declaration of the <!id>HeavyDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// <!id>HeavyDecayer<!!id> is a class that defines all the general routines 
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
#include <fstream>

namespace Herwig {

using namespace ThePEG;

class HeavyDecayer: public Decayer {

public:

  inline HeavyDecayer();
  inline HeavyDecayer(const HeavyDecayer &);
  virtual ~HeavyDecayer();
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

  static double VAWt(double*);//, double, double, double) const;
  // Weighting of phase space for V-A matrix elements

  //double PhaseSpaceWt() const;
  // Flat phase space weight (1.0)

  static ClassDescription<HeavyDecayer> initHeavyDecayer;
  static long lastAddedNumber;

  const HeavyDecayer & operator=(const HeavyDecayer &);
  //  Private and non-existent assignment operator.
};

}

namespace ThePEG {

template <>
struct BaseClassTrait<Herwig::HeavyDecayer,1> {
  typedef HandlerBase NthBase;
};

template <>
struct ClassTraits<Herwig::HeavyDecayer>: public ClassTraitsBase<Herwig::HeavyDecayer> {
  static string className() { return "/Herwig++/HeavyDecayer"; }
  static string library() { return "libHwDecay.so"; }
};

}

#include "HeavyDecayer.icc"

#endif /* HERWIG_HeavyDecayer_H */
