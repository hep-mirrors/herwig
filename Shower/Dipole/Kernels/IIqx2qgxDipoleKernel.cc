// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIqx2qgxDipoleKernel class.
//

#include "IIqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IIqx2qgxDipoleKernel::IIqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

IIqx2qgxDipoleKernel::~IIqx2qgxDipoleKernel() {}

IBPtr IIqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IIqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IIqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool IIqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    emitter(a) == sk.emitter(b) &&
    emission(a) == sk.emission(b) &&
    a.emitterPDF() == b.emitterPDF() &&
    a.spectatorData() == b.spectatorData() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr IIqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  return ind.emitterData();
}

tcPDPtr IIqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IIqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IIqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());
  double x = z*(1.-z)/(1.-z+ratio);

  ret *= (!strictLargeN() ? 4./3. : 3./2.) * ( (1.+sqr(x))/(1.-x) );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IIqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IIqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IIqx2qgxDipoleKernel> IIqx2qgxDipoleKernel::initIIqx2qgxDipoleKernel;
// Definition of the static class description member.

void IIqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<IIqx2qgxDipoleKernel> documentation
    ("IIqx2qgxDipoleKernel");

}

