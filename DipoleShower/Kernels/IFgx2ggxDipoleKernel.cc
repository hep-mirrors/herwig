// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFgx2ggxDipoleKernel class.
//

#include "IFgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFgx2ggxDipoleKernel::IFgx2ggxDipoleKernel() 
  : DipoleSplittingKernel() {}

IFgx2ggxDipoleKernel::~IFgx2ggxDipoleKernel() {}

IBPtr IFgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());
  double x = ( z*(1.-z) - ratio ) / ( 1. - z - ratio );
  double u = ratio / (1.-z);

  //double rho = (1.+x+u)/x;
  //ret *= 6. * ( 1./(1.-x+u) + (1.-x)/x - 1. + x*(1.-x) +u/rho );
  ret *= 3. * ( 1./(1.-x+u) + (1.-x)/x - 1. + x*(1.-x) );

  return ret;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFgx2ggxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFgx2ggxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFgx2ggxDipoleKernel> IFgx2ggxDipoleKernel::initIFgx2ggxDipoleKernel;
// Definition of the static class description member.

void IFgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<IFgx2ggxDipoleKernel> documentation
    ("IFgx2ggxDipoleKernel");

}

