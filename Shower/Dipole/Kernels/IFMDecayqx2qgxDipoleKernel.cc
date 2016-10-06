// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMDecayqx2qgxDipoleKernel class.
//

#include "IFMDecayqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMDecayqx2qgxDipoleKernel::IFMDecayqx2qgxDipoleKernel() : DipoleSplittingKernel() {}

IFMDecayqx2qgxDipoleKernel::~IFMDecayqx2qgxDipoleKernel() {}

IBPtr IFMDecayqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFMDecayqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFMDecayqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    ind.incomingDecayEmitter() && !ind.incomingDecaySpectator() &&
    // Only for top decay
    abs(ind.emitterData()->id()) == 6 &&
    // Matches to the kernel declared in a .in file for the given emitter flavour
    abs(ind.emitterData()->id()) == abs(flavour()->id()) &&
    !(ind.emitterData()->mass() == ZERO) && 
    // Initial state here refers to the entire event
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFMDecayqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
						     const DipoleSplittingKernel& sk,
						     const DipoleIndex& b) const {
  
  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    abs(sk.emitter(b)->id()) == 6 &&
    abs(sk.emitter(b)->mass()) == abs(emitter(a)->mass()) &&
    abs(sk.spectator(b)->mass()) == abs(spectator(a)->mass());

}

tcPDPtr IFMDecayqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {

  assert(flavour());
  assert(abs(flavour()->id()) == 6);

  return ind.emitterData()->id() > 0 ?
    (tcPDPtr) flavour() : (tcPDPtr) flavour()->CC();
}


tcPDPtr IFMDecayqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}


tcPDPtr IFMDecayqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFMDecayqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
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
  Energy2 mi2 = sqr(split.spectatorData()->mass());
  Energy2 ma2 = sqr(split.emitterData()->mass());
  
  // Recoil system mass
  Energy2 mk2 = sqr(split.recoilMass());

  Energy2 Qijk = sqr(split.scale());
  Energy2 sbar = Qijk - mi2 - mk2;

  // Calculate y
  double y = (1./(sbar*zPrime*(1.-zPrime)))*( sqr(pt) + sqr(1.-zPrime)*mi2 );

  double mui2 = mi2 / Qijk;
  double muk2 = mk2 / Qijk;
  double mua2 = ma2 / Qijk;
  double sbarMod = 1. - mui2 - muk2;
 
  if( sqr(2.*muk2+sbarMod*(1.-y))-4.*muk2 < 0. ){
    generator()->logWarning( Exception()
                            << "error in IFMDecayqx2qgxDipoleKernel::evaluate -- " <<
                            "muk2 " << muk2 << "  mui2 " << mui2 << "  y " << y
                            << Exception::warning );
    return 0.0;
  }

  double vijk = sqrt( sqr(2.*muk2+sbarMod*(1.-y))-4.*muk2 ) / (sbarMod*(1.-y));
  double vbar = sqrt( 1.+sqr(mui2)+sqr(muk2)-2.*(mui2+muk2+mui2*muk2) ) / sbarMod;


  ret *= (!strictLargeN() ? 4./3. : 3./2.)*( 2.*(2.*mui2 + 2.*y*sbarMod + sbarMod)/((1.+y)*sbarMod-z*(1.-y)*sbarMod) - (vbar/vijk)*(2. + 2.*mua2/((1.-z*(1.-y))*sbarMod)) );

  return ret > 0. ? ret : 0.;
  
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMDecayqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFMDecayqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFMDecayqx2qgxDipoleKernel> IFMDecayqx2qgxDipoleKernel::initIFMDecayqx2qgxDipoleKernel;
// Definition of the static class description member.

void IFMDecayqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<IFMDecayqx2qgxDipoleKernel> documentation
    ("IFMDecayqx2qgxDipoleKernel");

}

