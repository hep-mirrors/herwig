// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFqx2qgxDipoleKernel class.
//

#include "FFqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFqx2qgxDipoleKernel::FFqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

FFqx2qgxDipoleKernel::~FFqx2qgxDipoleKernel() {}

IBPtr FFqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() == ZERO &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}
#ifndef NDEBUG
bool FFqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
#else
bool FFqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& ,
#endif
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    abs(sk.emitter(b)->id()) < 6 &&
    sk.emitter(b)->mass() == ZERO;
       

}

tcPDPtr FFqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  return ind.emitterData();
}

tcPDPtr FFqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double y = sqr(split.lastPt() / split.scale()) / (z*(1.-z));

  ret *= (!strictLargeN() ? 4./3. : 3./2.)*( 2./(1.-z*(1.-y)) - (1.+z) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FFqx2qgxDipoleKernel::generatePhi(const DipoleSplittingInfo&, const RhoDMatrix&) const {
  
  // No dependence on the spin density matrix,
  // dependence on off-diagonal terms cancels.
  return {{ {0, 1.} }};
}

DecayMEPtr FFqx2qgxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {

  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt( 1./(1.-z) );
  double v_AP_ppm = -z/sqrt(1.-z);

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -1 or -1/2, 1=+1/2, 2 = +1
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,1,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(1,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = 0.;
  (*kernelPhiDep)(1,0,2) = 0.;
  (*kernelPhiDep)(0,1,2) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FFqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFqx2qgxDipoleKernel> FFqx2qgxDipoleKernel::initFFqx2qgxDipoleKernel;
// Definition of the static class description member.

void FFqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<FFqx2qgxDipoleKernel> documentation
    ("FFqx2qgxDipoleKernel");

}

