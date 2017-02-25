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
  useThisKernel() &&
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
  double mua2 = sqr( split.spectatorData()->mass() / split.scale() );
  double mui2 = sqr(split.emitterData()->mass() / split.scale());
  double mu2 = mui2;
  // Recoil system mass
  double muj2 = sqr(split.recoilMass() / split.scale());
  double bar = 1. - mui2 - mu2 - muj2;

  // Calculate y
  double y = (sqr(pt)/sqr(split.scale()) + sqr(1.-zPrime)*mui2 + sqr(zPrime)*mu2) / (bar*zPrime*(1.-zPrime));

  if( sqr(2.*muj2+bar*(1.-y))-4.*muj2 < 0. ){
    generator()->logWarning( Exception()
                            << "error in FIMDecaygx2qqxDipoleKernel::evaluate -- " <<
                            "muj2 " << muj2 << "  mui2 " << mui2 << "  y " << y << Exception::warning );

    return 0.0;
  }

  double vijk = sqrt( sqr(2.*muj2+bar*(1.-y))-4.*muj2 ) / (bar*(1.-y));
  double viji = sqrt( sqr(bar*y)-4.*sqr(mui2) ) / (bar*y+2.*mui2);
  double vbar = sqrt( 1.+sqr(mui2)+sqr(muj2)-2.*(mui2+muj2+mui2*muj2) ) / bar;
  
  double zp = 0.5*(1.+viji*vijk);
  double zm = 0.5*(1.-viji*vijk);

  // how to choose kappa?
  double kappa = 0.;

  ret *= 0.25 / vijk * ( 1. - 2.*( z*(1.-z) - (1.-kappa)*zp*zm - kappa*mui2/(2*mui2+bar*y) ) )
    + (!strictLargeN() ? 4./3. : 3./2.) * y/(1.-z*(1.-y)) * (-1.*(vbar/vijk)*2.*mua2/((1.-z*(1.-y))*bar) );
  
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
