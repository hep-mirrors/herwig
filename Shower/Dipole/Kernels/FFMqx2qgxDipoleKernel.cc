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
    abs(ind.emitterData()->id()) < 7  &&
    // 2012-05-01
    abs(ind.emitterData()->id()) == abs(flavour()->id()) &&
    !( ind.emitterData()->mass() == ZERO &&
       ind.spectatorData()->mass() == ZERO ) &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
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
  //  return ind.emitterData();
  //  cout << "q2qg kernel: flavour " << flavour()->id() << endl << flush;
  assert(flavour());
  assert(abs(flavour()->id())<7);
  //  cout << "  assertions passed" << endl << flush;
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
  
  // masses
  double muQ2 = sqr( split.emitterData()->mass() / split.scale() );
  double muj2 = sqr( split.spectatorData()->mass() / split.scale() );

  double z = split.lastZ();
  double y = ( sqr(split.lastPt()/split.scale()) + muQ2*sqr(1.-z) ) /
    (z*(1.-z)) / (1.-muQ2-muj2);

  double vijk = sqrt( sqr(2.*muj2+(1.-muQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-muQ2-muj2)*(1.-y));
  double vbar = sqrt( 1.+sqr(muQ2)+sqr(muj2)-2.*(muQ2+muj2+muQ2*muj2) ) / (1.-muQ2-muj2);

  ret *= (4./3.)*( 2./(1.-z*(1.-y)) - vbar/vijk*( 1.+z + muQ2*2./(y*(1.-muQ2-muj2)) ) );
  
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

