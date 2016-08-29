// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIqx2gqxDipoleKernel class.
//

#include "IIqx2gqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IIqx2gqxDipoleKernel::IIqx2gqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IIqx2gqxDipoleKernel::~IIqx2gqxDipoleKernel() {}

IBPtr IIqx2gqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IIqx2gqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IIqx2gqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool IIqx2gqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    a.emitterData() == b.emitterData() &&
    emitter(a) == sk.emitter(b) &&
    a.emitterPDF() == b.emitterPDF() &&
    a.spectatorData() == b.spectatorData() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr IIqx2gqxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IIqx2gqxDipoleKernel::emission(const DipoleIndex& ind) const {
  return ind.emitterData()->CC();
}

tcPDPtr IIqx2gqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IIqx2gqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());
  double x = z*(1.-z)/(1.-z+ratio);

  ret *= .5 * ( 1.-2.*x*(1.-x) );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IIqx2gqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IIqx2gqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IIqx2gqxDipoleKernel> IIqx2gqxDipoleKernel::initIIqx2gqxDipoleKernel;
// Definition of the static class description member.

void IIqx2gqxDipoleKernel::Init() {

  static ClassDocumentation<IIqx2gqxDipoleKernel> documentation
    ("IIqx2gqxDipoleKernel");

}

