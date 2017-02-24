// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIgx2ggxDipoleKernel class.
//

#include "IIgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IIgx2ggxDipoleKernel::IIgx2ggxDipoleKernel() 
  : DipoleSplittingKernel() {}

IIgx2ggxDipoleKernel::~IIgx2ggxDipoleKernel() {}

IBPtr IIgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IIgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IIgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool IIgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    a.emitterPDF() == b.emitterPDF() &&
    a.spectatorData() == b.spectatorData() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr IIgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IIgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IIgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IIgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());
  double x = z*(1.-z)/(1.-z+ratio);

  ret *= 3. * ( x/(1.-x) + (1.-x)/x +x*(1.-x) );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IIgx2ggxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IIgx2ggxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IIgx2ggxDipoleKernel> IIgx2ggxDipoleKernel::initIIgx2ggxDipoleKernel;
// Definition of the static class description member.

void IIgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<IIgx2ggxDipoleKernel> documentation
    ("IIgx2ggxDipoleKernel");

}

