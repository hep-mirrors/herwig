// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIqx2qgxDipoleKernel class.
//

#include "FIqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIqx2qgxDipoleKernel::FIqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

FIqx2qgxDipoleKernel::~FIqx2qgxDipoleKernel() {}

IBPtr FIqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6 &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    abs(sk.emitter(b)->id()) < 6 &&
    sk.emitter(b)->mass() == ZERO &&
    a.spectatorPDF() == b.spectatorPDF();
       

}


tcPDPtr FIqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  return ind.emitterData();
}

tcPDPtr FIqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double x = 1. / ( 1. + sqr(split.lastPt()/split.scale()) / (z*(1.-z)) );

  ret *= (!strictLargeN() ? 4./3. : 3./2.) * ( 2./(1.-z+(1.-x)) -(1.+z) + (1.-x)*(1.+3.*x*z) );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIqx2qgxDipoleKernel> FIqx2qgxDipoleKernel::initFIqx2qgxDipoleKernel;
// Definition of the static class description member.

void FIqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<FIqx2qgxDipoleKernel> documentation
    ("FIqx2qgxDipoleKernel");

}

