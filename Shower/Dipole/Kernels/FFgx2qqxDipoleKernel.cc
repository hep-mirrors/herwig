// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFgx2qqxDipoleKernel class.
//

#include "FFgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFgx2qqxDipoleKernel::FFgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

FFgx2qqxDipoleKernel::~FFgx2qqxDipoleKernel() {}

IBPtr FFgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() == ZERO &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool FFgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    sk.emitter(b)->mass() == ZERO;
       
}


tcPDPtr FFgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr FFgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour()->CC();
}

tcPDPtr FFgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();

  ret *= .25 * ( 1. - 2.*z*(1.-z) );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FFgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFgx2qqxDipoleKernel> FFgx2qqxDipoleKernel::initFFgx2qqxDipoleKernel;
// Definition of the static class description member.

void FFgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FFgx2qqxDipoleKernel> documentation
    ("FFgx2qqxDipoleKernel");

}

