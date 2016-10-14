// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFgx2qqxDipoleKernel class.
//

#include "IFgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFgx2qqxDipoleKernel::IFgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IFgx2qqxDipoleKernel::~IFgx2qqxDipoleKernel() {}

IBPtr IFgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() == ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    flavour() == sk.flavour() &&
    abs(flavour()->id()) < 6 &&
    flavour()->mass() == ZERO &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IFgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IFgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());

  double rho = 1. - 4.*ratio*z*(1.-z)/sqr(1.-z+ratio);
  double x = 0.5*((1.-z+ratio)/ratio)*(1.-sqrt(rho));

  ret *= 0.5 * (!strictLargeN() ? 4./3. : 3./2.) * ( x + 2.*(1.-x)/x );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFgx2qqxDipoleKernel> IFgx2qqxDipoleKernel::initIFgx2qqxDipoleKernel;
// Definition of the static class description member.

void IFgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<IFgx2qqxDipoleKernel> documentation
    ("IFgx2qqxDipoleKernel");

}

