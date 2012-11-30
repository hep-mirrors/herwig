// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMgx2qqxDipoleKernel class.
//

#include "FIMgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMgx2qqxDipoleKernel::FIMgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

FIMgx2qqxDipoleKernel::~FIMgx2qqxDipoleKernel() {}

IBPtr FIMgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() != ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIMgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    emitter(a)->mass() == sk.emitter(b)->mass() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr FIMgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() != ZERO);
  return flavour();
}

tcPDPtr FIMgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() != ZERO);
  return flavour()->CC();
}

tcPDPtr FIMgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

// TODO
// assure split.scale() is sqrt(sbar)
double FIMgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  // mi=m=mQ, Mi=0, mj=Mj=0
  Energy2 mQ2 = sqr(split.emitterData()->mass());

  double z = split.lastZ();
  double x = 1./ ( 1. +
		   ( sqr(split.lastPt()) + mQ2 ) /
		   ( z*(1.-z) * sqr(split.scale()) ) );

  double muQ2 = x * mQ2/sqr(split.scale());

  double zm = .5 * ( 1. - sqrt( 1. - 4.*sqr(muQ2/(1.-x)) ) );
  double zp = .5 * ( 1. - sqrt( 1. - 4.*sqr(muQ2/(1.-x)) ) );

  ret *= .25 * (1.-2.*(zp-z)*(z-zm));

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIMgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIMgx2qqxDipoleKernel> FIMgx2qqxDipoleKernel::initFIMgx2qqxDipoleKernel;
// Definition of the static class description member.

void FIMgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FIMgx2qqxDipoleKernel> documentation
    ("FIMgx2qqxDipoleKernel");

}

