// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMDecaygx2qqxDipoleKernel class.
//

#include "FIMDecaygx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMDecaygx2qqxDipoleKernel::FIMDecaygx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

FIMDecaygx2qqxDipoleKernel::~FIMDecaygx2qqxDipoleKernel() {}

IBPtr FIMDecaygx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMDecaygx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMDecaygx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    ind.incomingDecaySpectator() && !ind.incomingDecayEmitter() &&
    ind.emitterData()->id() == ParticleID::g &&
    !(ind.spectatorData()->mass() == ZERO) &&
    // Initial state here refers to the entire event
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool FIMDecaygx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
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


tcPDPtr FIMDecaygx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6);
  return flavour();
}

tcPDPtr FIMDecaygx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6);
  return flavour()->CC();
}

tcPDPtr FIMDecaygx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIMDecaygx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);

  // These are the physical variables as used in the standard form of the kernel (i.e. do not redefine variables or kernel)
  double z = split.lastZ();
  Energy pt = split.lastPt();

  // Need zPrime to calculate y, 
  // TODO: Should just store y in the dipole splitting info everywhere anyway!!!
  // The only value stored in dInfo.lastSplittingParameters() should be zPrime
  //assert(split.lastSplittingParameters().size() == 1 );
  double zPrime = split.lastSplittingParameters()[0];

  // Construct mass squared variables
  Energy2 mi2 = sqr(split.emitterData()->mass());

  // Recoil system mass
  Energy2 mk2 = sqr(split.recoilMass());

  Energy2 Qijk = sqr(split.scale());
  Energy2 sbar = Qijk - 2.*mi2 - mk2;

  // Calculate y
  double y = (1./(sbar*zPrime*(1.-zPrime)))*( sqr(pt) + (1. - 2.*zPrime + 2.*sqr(zPrime))*mi2 );

  // Scaled masses
  double mui2 = mi2 / Qijk;
  double muk2 = mk2 / Qijk;
  double sbarMod = 1. - 2.*mui2 - muk2;

  if( sqr(2.*muk2+sbarMod*(1.-y))-4.*muk2 < 0. ){
    cerr << "error in FIMDecaygx2qqxDipoleKernel::evaluate -- " <<
      "muk2 " << muk2 << "  mui2 " << mui2 << "  y " << y << endl;
    return 0.0;
  }

  double vijk = sqrt( sqr(2.*muk2+sbarMod*(1.-y))-4.*muk2 ) / (sbarMod*(1.-y));
  double viji = sqrt( sqr(sbarMod*y)-4.*sqr(mui2) ) / (sbarMod*y+2.*mui2);

  double zp = 0.5*(1.+viji*vijk);
  double zm = 0.5*(1.-viji*vijk);

  // how to choose kappa?
  double kappa = 0.;

  ret *= 0.25 / vijk *
    ( 1. - 2.*( z*(1.-z) - (1.-kappa)*zp*zm - kappa*mui2/(2*mui2+sbarMod*y) ) );
    
  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMDecaygx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIMDecaygx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIMDecaygx2qqxDipoleKernel> FIMDecaygx2qqxDipoleKernel::initFIMDecaygx2qqxDipoleKernel;
// Definition of the static class description member.

void FIMDecaygx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FIMDecaygx2qqxDipoleKernel> documentation
    ("FIMDecaygx2qqxDipoleKernel");

}
