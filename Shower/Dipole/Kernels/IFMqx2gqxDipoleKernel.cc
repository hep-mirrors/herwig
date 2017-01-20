// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMqx2gqxDipoleKernel class.
//

#include "IFMqx2gqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMqx2gqxDipoleKernel::IFMqx2gqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IFMqx2gqxDipoleKernel::~IFMqx2gqxDipoleKernel() {}

IBPtr IFMqx2gqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFMqx2gqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFMqx2gqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() != ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFMqx2gqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    a.emitterData() == b.emitterData() &&
    emitter(a) == sk.emitter(b) &&
    spectator(a)->mass() == sk.spectator(b)->mass() &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFMqx2gqxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFMqx2gqxDipoleKernel::emission(const DipoleIndex& ind) const {
  return ind.emitterData()->CC();
}

tcPDPtr IFMqx2gqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFMqx2gqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  Energy pt = split.lastPt();
  double ratio = sqr(pt/split.scale());
  double muk2 = sqr(split.spectatorData()->mass()/split.scale());
  
// Calculate x
    double rho = 1. - 4.*ratio*(1.-muk2)*z*(1.-z)/sqr(1.-z+ratio);
    double x = 0.5*((1.-z+ratio)/(ratio*(1.-muk2))) * (1. - sqrt(rho));


  ret *= .5 * ( 1.-2.*x*(1.-x)  );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMqx2gqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFMqx2gqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFMqx2gqxDipoleKernel> IFMqx2gqxDipoleKernel::initIFMqx2gqxDipoleKernel;
// Definition of the static class description member.

void IFMqx2gqxDipoleKernel::Init() {

  static ClassDocumentation<IFMqx2gqxDipoleKernel> documentation
    ("IFMqx2gqxDipoleKernel");

}

