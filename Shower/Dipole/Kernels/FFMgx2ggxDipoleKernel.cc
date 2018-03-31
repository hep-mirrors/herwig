// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMgx2ggxDipoleKernel class.
//

#include "FFMgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFMgx2ggxDipoleKernel::FFMgx2ggxDipoleKernel() 
  : DipoleSplittingKernel(){}

FFMgx2ggxDipoleKernel::~FFMgx2ggxDipoleKernel() {}

IBPtr FFMgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFMgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFMgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() != ZERO &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator() &&
    !ind.incomingDecayEmitter() && !ind.incomingDecaySpectator();

}

bool FFMgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {
    
  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    abs(spectator(a)->mass()) == abs(sk.spectator(b)->mass());

}

tcPDPtr FFMgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFMgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFMgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFMgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);

  // These are the physical variables as used in the 
  // standard form of the kernel (i.e. do not redefine variables or kernel)
  const double z = split.lastZ();
  const Energy pt = split.lastPt();

  // Need zPrime to calculate y, 
  // TODO: Should just store y in the dipole splitting info everywhere anyway!!!
  // The only value stored in dInfo.lastSplittingParameters() should be zPrime

  const double zPrime = split.lastSplittingParameters()[0];
  
  // Construct mass squared variables

  // Construct mass squared variables
  double muj2 = sqr(split.spectatorMass() / split.scale());
  double bar = 1. - muj2;

  // Calculate y
  double y = sqr(pt)/sqr(split.scale()) / (bar*zPrime*(1.-zPrime));

  double vijk = sqrt( sqr(2.*muj2+bar*(1.-y))-4.*muj2 ) / (bar*(1.-y));
  double viji = 1.;

  double zp = 0.5*(1.+viji*vijk);
  double zm = 0.5*(1.-viji*vijk);

  // how to choose kappa?
  double kappa = 0.;

  double S1=1./(1.-z*(1.-y));
  double S2=1./(1.-(1.-z)*(1.-y));
  double NS=(z*(1.-z)-(1.-kappa)*zp*zm-2.)/vijk;
  
  if( theAsymmetryOption == 0 ){
    ret *= 3.*( S1 + 0.5 * NS);
  }else if ( theAsymmetryOption == 1 ){
    ret *= 3.*z*( S1 +S2 + NS );
  }else{
    ret *= 3.*0.5*( S1 + S2 + NS );
  }

  return ret > 0. ? ret : 0.;  
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMgx2ggxDipoleKernel::persistentOutput(PersistentOStream & os) const {
  os << theAsymmetryOption;
}

void FFMgx2ggxDipoleKernel::persistentInput(PersistentIStream & is, int) {
  is >> theAsymmetryOption;
}

ClassDescription<FFMgx2ggxDipoleKernel> FFMgx2ggxDipoleKernel::initFFMgx2ggxDipoleKernel;
// Definition of the static class description member.

void FFMgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<FFMgx2ggxDipoleKernel> documentation
    ("FFMgx2ggxDipoleKernel");
  
  static Parameter<FFMgx2ggxDipoleKernel,int> interfacetheAsymmetryOption
  ("AsymmetryOption",
   "The asymmetry option for final state gluon spliitings.",
   &FFMgx2ggxDipoleKernel::theAsymmetryOption, 1, 0, 0,
   false, false, Interface::lowerlim);

}

