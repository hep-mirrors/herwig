// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMqx2qgxDipoleKernel class.
//

#include "IFMqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMqx2qgxDipoleKernel::IFMqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

IFMqx2qgxDipoleKernel::~IFMqx2qgxDipoleKernel() {}

IBPtr IFMqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFMqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFMqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() != ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFMqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    emitter(a) == sk.emitter(b) &&
    emission(a) == sk.emission(b) &&
    spectator(a)->mass() == sk.spectator(b)->mass() &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFMqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  return ind.emitterData();
}

tcPDPtr IFMqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFMqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFMqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  Energy pt = split.lastPt();
  double ratio = sqr(pt/split.scale());
  double muk2 = sqr(split.spectatorData()->mass()/split.scale());

// Calculate x and u
    double rho = 1. - 4.*ratio*(1.-muk2)*z*(1.-z)/sqr(1.-z+ratio);
    double x = 0.5*((1.-z+ratio)/(ratio*(1.-muk2))) * (1. - sqrt(rho));
    double u = x*ratio / (1.-z);
  
// 19/01/2017 - SW: Removed finite term as its effect on
// the ratio of -ve/+ve kernels is small
    ret *= (!strictLargeN() ? 4./3. : 3./2.) * ( 2./(1.-x+u) - (1.+x) );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFMqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFMqx2qgxDipoleKernel> IFMqx2qgxDipoleKernel::initIFMqx2qgxDipoleKernel;
// Definition of the static class description member.

void IFMqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<IFMqx2qgxDipoleKernel> documentation
    ("IFMqx2qgxDipoleKernel");

}

