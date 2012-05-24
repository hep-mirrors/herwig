// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMgx2ggxDipoleKernel class.
//

#include "FFMgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFMgx2ggxDipoleKernel::FFMgx2ggxDipoleKernel() 
  : DipoleSplittingKernel() {}

FFMgx2ggxDipoleKernel::~FFMgx2ggxDipoleKernel() {}

IBPtr FFMgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFMgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFMgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    ind.emitterData()->id() == ParticleID::g &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator() &&
    // massless quarks shall be treated by massless FFgx2ggxDipoleKernel
    !( abs(ind.spectatorData()->id()) < 6 && ind.spectatorData()->mass() == ZERO );
}

bool FFMgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {
    
  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    abs(spectator(a)->id()) == abs(sk.spectator(b)->id());

}

tcPDPtr FFMgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFMgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFMgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFMgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);
  
  // masses
  double muj2 = sqr( split.spectatorData()->mass() / split.scale() );
  
  double z = split.lastZ();
  double y = sqr(split.lastPt() / split.scale()) / (z*(1.-z)) / (1.-muj2);

  double vijk = sqrt( sqr(2.*muj2+(1.-muj2)*(1.-y))-4.*muj2 ) / ((1.-muj2)*(1.-y));
  double viji = 1.;

  double zp = 0.5*(1.+viji*vijk);
  double zm = 0.5*(1.-viji*vijk);

  // free parameter of the subtraction scheme
  double kappa = 0.;

  ret *= 3.*(1./(1.-z*(1.-y))+1./(1.-(1.-z)*(1.-y)) + (z*(1.-z)-(1.-kappa)*zp*zm-2.)/vijk);
  
  return ret > 0. ? ret : 0.;
  
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMgx2ggxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FFMgx2ggxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFMgx2ggxDipoleKernel> FFMgx2ggxDipoleKernel::initFFMgx2ggxDipoleKernel;
// Definition of the static class description member.

void FFMgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<FFMgx2ggxDipoleKernel> documentation
    ("FFMgx2ggxDipoleKernel");

}

