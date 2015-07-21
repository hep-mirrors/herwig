// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIgx2qqxDipoleKernel class.
//

#include "IIgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IIgx2qqxDipoleKernel::IIgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IIgx2qqxDipoleKernel::~IIgx2qqxDipoleKernel() {}

IBPtr IIgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IIgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IIgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() == ZERO &&
    ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool IIgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    flavour() == sk.flavour() &&
    abs(flavour()->id()) < 6 &&
    flavour()->mass() == ZERO &&
    a.emitterPDF() == b.emitterPDF() &&
    a.spectatorData() == b.spectatorData() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr IIgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IIgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IIgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IIgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());
  double x = z*(1.-z)/(1.-z+ratio);

  ret *= 0.5 * (!strictLargeN() ? 4./3. : 3./2.) * ( 1./x +sqr(1.-x)/x );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IIgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IIgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IIgx2qqxDipoleKernel> IIgx2qqxDipoleKernel::initIIgx2qqxDipoleKernel;
// Definition of the static class description member.

void IIgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<IIgx2qqxDipoleKernel> documentation
    ("IIgx2qqxDipoleKernel");

}

