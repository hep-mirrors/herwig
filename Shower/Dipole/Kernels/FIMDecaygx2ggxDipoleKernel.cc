// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMDecaygx2ggxDipoleKernel class.
//

#include "FIMDecaygx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMDecaygx2ggxDipoleKernel::FIMDecaygx2ggxDipoleKernel() 
  : DipoleSplittingKernel(),theSymmetryFactor(0.5){}

FIMDecaygx2ggxDipoleKernel::~FIMDecaygx2ggxDipoleKernel() {}

IBPtr FIMDecaygx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMDecaygx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMDecaygx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    ind.incomingDecaySpectator() && !ind.incomingDecayEmitter() && 
    ind.emitterData()->id() == ParticleID::g &&
    !(ind.spectatorData()->mass() == ZERO) &&
    // Initial state here refers to the entire event
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool FIMDecaygx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    sk.emitter(b)->id() == ParticleID::g &&
    abs(sk.spectator(b)->mass()) == abs(spectator(a)->mass());

}


tcPDPtr FIMDecaygx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIMDecaygx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIMDecaygx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIMDecaygx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

    double ret = alphaPDF(split);

  // These are the physical variables as used in the standard form of the kernel (i.e. do not redefine variables or kernel)
  double z = split.lastZ();
  Energy pt = split.lastPt();

  // Need zPrime to calculate y, 
  // TODO: Should just store y in the dipole splitting info everywhere anyway!!!
  // The only value stored in dInfo.lastSplittingParameters() should be zPrime
  //assert(split.lastSplittingParameters().size() == 1 );
  double zPrime = split.lastSplittingParameters()[0];

  // Recoil system mass
  Energy2 mk2 = sqr(split.recoilMass());

  Energy2 Qijk = sqr(split.scale());
  Energy2 sbar = Qijk - mk2;

  // Calculate y
  double y = (1./(sbar*zPrime*(1.-zPrime)))*( sqr(pt) );

  // Scaled masses
  double muk2 = mk2 / Qijk;
  double sbarMod = 1. - muk2;


  if( sqr(2.*muk2+sbarMod*(1.-y))-4.*muk2 < 0. ){
    generator()->logWarning( Exception()
                            << "error in FIMDecaygx2ggxDipoleKernel::evaluate -- " <<
                            "muk2 " << muk2 << "  y " << y << Exception::warning );
    return 0.0;
  }

  double vijk = sqrt( sqr(2.*muk2+sbarMod*(1.-y))-4.*muk2 ) / (sbarMod*(1.-y));
  double viji = 1.;

  double zp = 0.5*(1.+viji*vijk);
  double zm = 0.5*(1.-viji*vijk);

  // how to choose kappa?
  double kappa = 0.;

  ret *= theSymmetryFactor*3.*( (2.*y + 1.)/((1.+y)-z*(1.-y)) + (2.*y + 1.)/((1.+y)-(1.-z)*(1.-y)) + (1./vijk)*( z*(1.-z) - (1.-kappa)*zp*zm - 2. ) ) ;

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMDecaygx2ggxDipoleKernel::persistentOutput(PersistentOStream & os) const {

  os<<theSymmetryFactor;
}

void FIMDecaygx2ggxDipoleKernel::persistentInput(PersistentIStream & is, int) {

  is>>theSymmetryFactor;
}

ClassDescription<FIMDecaygx2ggxDipoleKernel> FIMDecaygx2ggxDipoleKernel::initFIMDecaygx2ggxDipoleKernel;
// Definition of the static class description member.

void FIMDecaygx2ggxDipoleKernel::Init() {

  static ClassDocumentation<FIMDecaygx2ggxDipoleKernel> documentation
    ("FIMDecaygx2ggxDipoleKernel");

  static Parameter<FIMDecaygx2ggxDipoleKernel,double> interfaceSymmetryFactor
    ("SymmetryFactor",
     "The symmetry factor for final state gluon splittings.",
     &FIMDecaygx2ggxDipoleKernel::theSymmetryFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);
}
