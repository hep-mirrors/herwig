// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIgx2qqxDipoleKernel class.
//

#include "FIgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIgx2qqxDipoleKernel::FIgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

FIgx2qqxDipoleKernel::~FIgx2qqxDipoleKernel() {}

IBPtr FIgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() == ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    // sk.emitter(b)->mass() == ZERO &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr FIgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr FIgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour()->CC();
}

tcPDPtr FIgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();

  ret *= .25 * (1.-2.*z*(1.-z));

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIgx2qqxDipoleKernel> FIgx2qqxDipoleKernel::initFIgx2qqxDipoleKernel;
// Definition of the static class description member.

void FIgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FIgx2qqxDipoleKernel> documentation
    ("FIgx2qqxDipoleKernel");

}

