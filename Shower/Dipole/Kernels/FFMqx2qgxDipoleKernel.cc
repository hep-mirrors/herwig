// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMqx2qgxDipoleKernel class.
//

#include "FFMqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFMqx2qgxDipoleKernel::FFMqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

FFMqx2qgxDipoleKernel::~FFMqx2qgxDipoleKernel() {}

IBPtr FFMqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFMqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFMqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 7  &&
    // 2012-05-01
    abs(ind.emitterData()->id()) == abs(flavour()->id()) &&
    !( ind.emitterData()->mass() == ZERO &&
       ind.spectatorData()->mass() == ZERO ) &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator() &&
    !ind.incomingDecayEmitter() && !ind.incomingDecaySpectator();
}

bool FFMqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {
  
  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    abs(sk.emitter(b)->id()) < 7 &&
    abs(sk.emitter(b)->mass()) == abs(emitter(a)->mass()) &&
    abs(sk.spectator(b)->mass()) == abs(spectator(a)->mass());

}

// 2012-05-01
tcPDPtr FFMqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  assert(flavour());
  assert(abs(flavour()->id())<7);
  return ind.emitterData()->id() > 0 ?
    (tcPDPtr) flavour() : (tcPDPtr) flavour()->CC();
}

tcPDPtr FFMqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFMqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFMqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);

  // These are the physical variables as used in the standard 
  // form of the kernel (i.e. do not redefine variables or kernel)
  double z = split.lastZ();
  Energy pt = split.lastPt();

  // Need zPrime to calculate y, 
  // TODO: Should just store y in the dipole splitting info everywhere anyway!!!
  // The only value stored in dInfo.lastSplittingParameters() should be zPrime
  //assert(split.lastSplittingParameters().size() == 1 );
  double zPrime = split.lastSplittingParameters()[0];

  // Construct mass squared variables
  double mui2 = sqr(split.emitterData()->mass() / split.scale());
  double muj2 = sqr(split.spectatorData()->mass() / split.scale());
  double bar = 1. - mui2 - muj2;

  // Calculate y
  double y = (sqr(pt)/sqr(split.scale()) + sqr(1.-zPrime)*mui2) / (bar*zPrime*(1.-zPrime));

  double vijk = sqrt( sqr(2.*muj2 + bar*(1.-y))-4.*muj2 ) / (bar*(1.-y));
  double vbar = sqrt( 1.+sqr(mui2)+sqr(muj2)-2.*(mui2+muj2+mui2*muj2) ) / bar;

  ret *= (!strictLargeN() ? 4./3. : 3./2.)*( 2./(1.-z*(1.-y)) - vbar/vijk*( 1.+z + mui2*2./(y*(1.-mui2-muj2)) ) );
  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FFMqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFMqx2qgxDipoleKernel> FFMqx2qgxDipoleKernel::initFFMqx2qgxDipoleKernel;
// Definition of the static class description member.

void FFMqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<FFMqx2qgxDipoleKernel> documentation
    ("FFMqx2qgxDipoleKernel");

}

