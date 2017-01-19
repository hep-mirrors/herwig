// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMgx2qqxDipoleKernel class.
//

#include "IFMgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMgx2qqxDipoleKernel::IFMgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IFMgx2qqxDipoleKernel::~IFMgx2qqxDipoleKernel() {}

IBPtr IFMgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFMgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFMgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() != ZERO &&
    flavour()->mass() == ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFMgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    flavour() == sk.flavour() &&
    abs(flavour()->id()) < 6 &&
    flavour()->mass() == ZERO &&
    spectator(a)->mass() == sk.spectator(b)->mass() &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFMgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IFMgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IFMgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFMgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  Energy pt = split.lastPt();
  double ratio = sqr(pt/split.scale());
  double muk2 = sqr(split.spectatorData()->mass()/split.scale());
  
// Calculate x and u
  double rho = 1. - 4.*ratio*(1.-muk2)*z*(1.-z)/sqr(1.-z+ratio);
  double x = 0.5*((1.-z+ratio)/(ratio*(1.-muk2))) * (1. - sqrt(rho));
  double u = x*ratio / (1.-z);

  // NOTE - The definition of muk used in the kinematics differs from that in CS
    double muk2CS = x*muk2;
    ret *= 0.5 * (!strictLargeN() ? 4./3. : 3./2.) *
      ( x + 2.*(1.-x)/x - 2.*muk2CS/x*u/(1.-u) );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFMgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFMgx2qqxDipoleKernel> IFMgx2qqxDipoleKernel::initIFMgx2qqxDipoleKernel;
// Definition of the static class description member.

void IFMgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<IFMgx2qqxDipoleKernel> documentation
    ("IFMgx2qqxDipoleKernel");

}

