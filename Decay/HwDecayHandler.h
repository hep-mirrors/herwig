// -*- C++ -*-
#ifndef HERWIG_HwDecayHandler_H
#define HERWIG_HwDecayHandler_H
//
// This is the declaration of the <!id>HwDecayHandler<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// The <!id>HwDecayHandler<!!id> is the base class of all handlers
// implementing the administration of decays of unstable particles. It
// is derived from the more general <!class>StepHandler<!!class>
// class, and overrides the <!id>handle<!!id> method, simply decays
// all unstable particle in the current step.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:StepHandler.html">StepHandler.h</a>
// <a href="http:CollisionHandler.html">CollisionHandler.h</a>
// <a href="http:SubProcessHandler.html">SubProcessHandler.h</a>
// 

#include "ThePEG/Handlers/DecayHandler.h"

using namespace ThePEG;

namespace Herwig {

class HwDecayHandler: public DecayHandler {

public:

  inline HwDecayHandler() : DecayHandler() {}
  inline HwDecayHandler(const HwDecayHandler &x) : DecayHandler(x) {}
  virtual ~HwDecayHandler();
  // Standard ctors and dtor

public:

  virtual void handle(PartialCollisionHandler & ch, const tPVector & tagged,
		      const Hint & hint)
    ThePEG_THROW_SPEC((Veto, Stop, Exception));
  // Look through all tagged particled and decay all unstable ones.

  static void Init();
  // Standard Init function
  void persistentOutput(PersistentOStream&) const;
  void persistentInput(PersistentIStream&, int);
protected:
  virtual IBPtr clone() const;
  virtual IBPtr fullclone() const;

  virtual void doupdate() throw(UpdateException);
  virtual void doinit() throw(InitException);
  virtual void dofinish();
  // Standard Interfaced virtual functions.

  virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.
  virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:
  static ClassDescription<HwDecayHandler> initHwDecayHandler;

  HwDecayHandler & operator=(const HwDecayHandler &);
  //  Private and non-existent assignment operator.

};
}

namespace ThePEG {
template <>
struct BaseClassTrait<Herwig::HwDecayHandler,1> {
  typedef DecayHandler NthBase;
};

template <>
struct ClassTraits<Herwig::HwDecayHandler>: public ClassTraitsBase<Herwig::HwDecayHandler> {
  static string className() { return "/Herwig++/HwDecayHandler"; }
};

}

#endif /* HERWIG_HwDecayHandler_H */
