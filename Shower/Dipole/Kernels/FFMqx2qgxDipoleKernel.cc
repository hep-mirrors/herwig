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
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 7  &&
    // 2012-05-01
    abs(ind.emitterData()->id()) == abs(flavour()->id()) &&
    !( ind.emitterData()->mass() == ZERO &&
       ind.spectatorData()->mass() == ZERO ) &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator() &&
    !ind.incomingDecayEmitter() && !ind.incomingDecaySpectator();
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
  assert(flavour());
  assert(abs(flavour()->id())<7);
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

  // Sudakov parameterisation variables,
  // needed to calculate y.
  double z = split.lastZ();
  Energy pt = split.lastPt();

  // Construct mass squared variables
  // Note for q->qg can use the emitterMass
  // (i.e. mass of emitter before splitting = mass of emitter after)
  Energy2 Qijk = sqr(split.scale());
  Energy2 mi2 = sqr(split.emitterMass());
  Energy2 mk2 = sqr(split.spectatorMass());
  Energy2 sbar = Qijk - mi2 - mk2;

  // Calculate y
  double y = (sqr(pt) + sqr(1.-z)*mi2) / sbar / z / (1.-z);

  // zi, used in dipole splitting kernel
  double zi = split.lastSplittingParameters()[0];

  double vijk = sqrt( sqr(2.*mk2 + sbar*(1.-y)) - 4.*mk2*Qijk ) / sbar / (1.-y);
  double vtilde = sqrt( sqr(Qijk) + sqr(mi2) + sqr(mk2)
                        - 2.*(mi2*Qijk + mk2*Qijk + mi2*mk2) ) / sbar;
  
  ret *= (!strictLargeN() ? 4./3. : 3./2.)*
    ( 2./(1.-zi*(1.-y)) - vtilde/vijk*( 1. + zi + 2.*mi2/sbar/y ) );
  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FFMqx2qgxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {
  
  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr FFMqx2qgxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {

  // Need variables for the AP kernels
  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();
  Energy mi = dInfo.emitterMass();

  // Altarelli-Parisi spin-indexed kernels:
  Energy den = sqrt(sqr(mi)*sqr(1.-z) + sqr(pt));
  double v_AP_ppp = pt / den / sqrt(1.-z);
  double v_AP_ppm = - z * v_AP_ppp ;
  double v_AP_pmp = mi*(1.-z)*sqrt(1.-z) / den ;

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = v_AP_pmp;

  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -1 or -1/2, 1=+1/2, 2 = +1
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,1,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(1,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm;
  (*kernelPhiDep)(1,0,2) = v_AP_pmp;
  (*kernelPhiDep)(0,1,2) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;

  return kernelPhiDep;
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

