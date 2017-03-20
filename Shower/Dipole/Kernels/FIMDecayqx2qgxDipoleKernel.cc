// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMDecayqx2qgxDipoleKernel class.
//

#include "FIMDecayqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMDecayqx2qgxDipoleKernel::FIMDecayqx2qgxDipoleKernel() : DipoleSplittingKernel() {}

FIMDecayqx2qgxDipoleKernel::~FIMDecayqx2qgxDipoleKernel() {}

IBPtr FIMDecayqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMDecayqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMDecayqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.incomingDecaySpectator() && !ind.incomingDecayEmitter() &&
    abs(ind.emitterData()->id()) < 7  &&
    // This line matches to the kernel declared in a .in file for the given emitter flavour
    abs(ind.emitterData()->id()) == abs(flavour()->id()) &&
    !(ind.spectatorData()->mass() == ZERO) && 
    // Initial state here refers to the entire event
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool FIMDecayqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
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

tcPDPtr FIMDecayqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {

  assert(flavour());
  assert(abs(flavour()->id()) < 7);

  return ind.emitterData()->id() > 0 ?
    (tcPDPtr) flavour() : (tcPDPtr) flavour()->CC();
}


tcPDPtr FIMDecayqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}


tcPDPtr FIMDecayqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}


double FIMDecayqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
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
  double mui2 = sqr(split.emitterData()->mass() / split.scale());
  // Recoil system mass
  double muj2 = sqr(split.recoilMass() / split.scale());
  double mua2 = sqr( split.spectatorData()->mass() / split.scale() );
  double bar = 1. - mui2 - muj2;
  
  // Calculate y
  double y = (sqr(pt)/sqr(split.scale()) + sqr(1.-zPrime)*mui2) / (bar*zPrime*(1.-zPrime));

  if( sqr(2.*muj2+bar*(1.-y))-4.*muj2 < 0. ){
    generator()->logWarning( Exception()
    << "error in FIMDecayqx2qgxDipoleKernel::evaluate -- " <<
    "muj2 " << muj2 << "  mui2 " << mui2 << "  y " << y << Exception::warning );
    return 0.0;
  }

  double vijk = sqrt( sqr(2.*muj2 + bar*(1.-y))-4.*muj2 ) / (bar*(1.-y));
  double vbar = sqrt( 1.+sqr(mui2)+sqr(muj2)-2.*(mui2+muj2+mui2*muj2) ) / bar;

  ret *=  (!strictLargeN() ? 4./3. : 3./2.) * ( ( 2.*(2.*mui2/bar + 2.*y + 1.)/((1.+y)-z*(1.-y)) - (vbar/vijk)*((1.+z) + 2.*mui2/(y*bar)) ) + y/(1.-z*(1.-y)) * ( 2.*(2.*mui2/bar + 2.*y + 1.)/((1.+y)-z*(1.-y)) - (vbar/vijk)*(2. + 2.*mua2/((1.-z*(1.-y))*bar)) ) );
  
  return ret > 0. ? ret : 0.;
  
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMDecayqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIMDecayqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIMDecayqx2qgxDipoleKernel> FIMDecayqx2qgxDipoleKernel::initFIMDecayqx2qgxDipoleKernel;
// Definition of the static class description member.

void FIMDecayqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<FIMDecayqx2qgxDipoleKernel> documentation
    ("FIMDecayqx2qgxDipoleKernel");

}

