// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFgx2qqxDipoleKernel class.
//

#include "FFgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFgx2qqxDipoleKernel::FFgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

FFgx2qqxDipoleKernel::~FFgx2qqxDipoleKernel() {}

IBPtr FFgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() == ZERO &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

#ifndef NDEBUG
bool FFgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
#else
bool FFgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& ,
#endif
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    sk.emitter(b)->mass() == ZERO;
       
}


tcPDPtr FFgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr FFgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour()->CC();
}

tcPDPtr FFgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();

  ret *= .25 * ( 1. - 2.*z*(1.-z) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FFgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppm = z;
  double v_AP_pmp = -(1.-z);

  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppm) + sqr(v_AP_pmp)) + 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);
  
  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*( sqr(v_AP_ppm) + sqr(v_AP_pmp) )/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );

  return distPhiDep;
}

DecayMEPtr FFgx2qqxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {
  
  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppm = z;
  double v_AP_pmp = -(1.-z);

  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());
  
  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = 0.;
  (*kernelPhiDep)(2,1,1) = 0.;
  (*kernelPhiDep)(0,0,1) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,1,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,1,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,1,1) = 0.;
  (*kernelPhiDep)(2,0,0) = 0.;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FFgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFgx2qqxDipoleKernel> FFgx2qqxDipoleKernel::initFFgx2qqxDipoleKernel;
// Definition of the static class description member.

void FFgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FFgx2qqxDipoleKernel> documentation
    ("FFgx2qqxDipoleKernel");

}

