// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMqx2qgxDipoleKernel class.
//

#include "FIMqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMqx2qgxDipoleKernel::FIMqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

FIMqx2qgxDipoleKernel::~FIMqx2qgxDipoleKernel() {}

IBPtr FIMqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 7 &&
    abs(ind.emitterData()->id())==abs(flavour()->id()) &&
    ind.emitterData()->mass() != ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIMqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    abs(sk.emitter(b)->id()) < 7 &&
    sk.emitter(b)->mass() == emitter(a)->mass() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr FIMqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  assert(flavour());
  assert(abs(flavour()->id())<7 && flavour()->mass() != ZERO);
  return ind.emitterData()->id() > 0 ?
    (tcPDPtr) flavour() : (tcPDPtr) flavour()->CC();
}

tcPDPtr FIMqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIMqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

// TODO
// split.scale() should be sqrt(sbar) = sqrt( Mi2 - Q2 ) !!!
double FIMqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  // Mi=mi=mQ, m=0, Mj=mj=0
  Energy2 mQ2 = sqr(split.emitterData()->mass());

  double z = split.lastZ();
  double x = 1. / ( 1. + 
		    ( sqr(split.lastPt()) + sqr(1.-z)*mQ2 ) /
		    ( z*(1.-z) * sqr(split.scale()) ) );

  // Simon has extra terms
  ret *= (!strictLargeN() ? 4./3. : 3./2.) *
    ( 2./(1.-z+(1.-x)) -(1.+z) - mQ2/sqr(split.scale()) * 2.*x/(1.-x) );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIMqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIMqx2qgxDipoleKernel> FIMqx2qgxDipoleKernel::initFIMqx2qgxDipoleKernel;
// Definition of the static class description member.

void FIMqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<FIMqx2qgxDipoleKernel> documentation
    ("FIMqx2qgxDipoleKernel");

}

