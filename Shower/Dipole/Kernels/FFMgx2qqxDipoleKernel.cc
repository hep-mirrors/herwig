// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMgx2qqxDipoleKernel class.
//

#include "FFMgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFMgx2qqxDipoleKernel::FFMgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

FFMgx2qqxDipoleKernel::~FFMgx2qqxDipoleKernel() {}

IBPtr FFMgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFMgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFMgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    ind.emitterData()->id() == ParticleID::g &&
    !( ind.spectatorData()->mass() == ZERO &&
       flavour()->mass() == ZERO ) &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool FFMgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {
  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    emitter(a)->id() == sk.emitter(b)->id() &&
    abs(sk.spectator(b)->mass()) == abs(spectator(a)->mass());
}

tcPDPtr FFMgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6);
  return flavour();
}

tcPDPtr FFMgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6);
  return flavour()->CC();
}

tcPDPtr FFMgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

/*
 * TODO: remove unnecessary if statement
 */
double FFMgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);
  
  // masses
  double muQ2 = sqr( split.emitterData()->mass() / split.scale() );
  double muj2 = sqr( split.spectatorData()->mass() / split.scale() );

  double z = split.lastZ();
  double y = ( sqr(split.lastPt()/split.scale()) + muQ2*(1.-2.*z+2.*z*z) ) /
    (z*(1.-z)) / (1.-2.*muQ2-muj2);
  
  // new: 2011-08-31
  // 2011-11-06: so far never happened
  if( sqr(2.*muj2+(1.-2.*muQ2-muj2)*(1.-y))-4.*muj2 < 0. ){
    cout << "error in FFMgx2qqxDipoleKernel::evaluate1 -- " <<
      "muj2 " << muj2 << "  muQ2 " << muQ2 << "  y " << y << endl;
    return 0.0;
  }
  
  double vijk = sqrt( sqr(2.*muj2+(1.-2.*muQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-2.*muQ2-muj2)*(1.-y));
  double viji = sqrt( sqr((1.-2.*muQ2-muj2)*y)-4.*sqr(muQ2) ) / ((1.-2.*muQ2-muj2)*y+2.*muQ2);

  double zp = 0.5*(1.+viji*vijk);
  double zm = 0.5*(1.-viji*vijk);

  // how to choose kappa??
  double kappa = 0.;

  ret *= 0.25 / vijk *
    ( 1. - 2.*( z*(1.-z) - (1.-kappa)*zp*zm - kappa*muQ2/(2*muQ2+(1.-2*muQ2-muj2)*y) ) );
    
  return ret > 0. ? ret : 0.;
  
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FFMgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFMgx2qqxDipoleKernel> FFMgx2qqxDipoleKernel::initFFMgx2qqxDipoleKernel;
// Definition of the static class description member.

void FFMgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FFMgx2qqxDipoleKernel> documentation
    ("FFMgx2qqxDipoleKernel");

}

