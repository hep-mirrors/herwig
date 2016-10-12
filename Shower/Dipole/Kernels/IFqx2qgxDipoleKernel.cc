// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFqx2qgxDipoleKernel class.
//

#include "IFqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFqx2qgxDipoleKernel::IFqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

IFqx2qgxDipoleKernel::~IFqx2qgxDipoleKernel() {}

IBPtr IFqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    emitter(a) == sk.emitter(b) &&
    emission(a) == sk.emission(b) &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  return ind.emitterData();
}

tcPDPtr IFqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());

  double rho = 1. - 4.*ratio*z*(1.-z)/sqr(1.-z+ratio);
  double x = 0.5*((1.-z+ratio)/ratio)*(1.-sqrt(rho));
  double u = 0.5*((1.-z+ratio)/(1.-z))*(1.-sqrt(rho));

  ret *= (!strictLargeN() ? 4./3. : 3./2.) * ( 2./(1.-x+u) - (1.+x) + u*(1.+3.*x*(1.-u) ) );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFqx2qgxDipoleKernel> IFqx2qgxDipoleKernel::initIFqx2qgxDipoleKernel;
// Definition of the static class description member.

void IFqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<IFqx2qgxDipoleKernel> documentation
    ("IFqx2qgxDipoleKernel");

}

