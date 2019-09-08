// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIgx2qqxDipoleKernel class.
//

#include "FIgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIgx2qqxDipoleKernel::FIgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

FIgx2qqxDipoleKernel::~FIgx2qqxDipoleKernel() {}

IBPtr FIgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() == ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    // sk.emitter(b)->mass() == ZERO &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr FIgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr FIgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour()->CC();
}

tcPDPtr FIgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();

  ret *= .25 * (1.-2.*z*(1.-z));

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FIgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

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

DecayMEPtr FIgx2qqxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {
  
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


void FIgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIgx2qqxDipoleKernel> FIgx2qqxDipoleKernel::initFIgx2qqxDipoleKernel;
// Definition of the static class description member.

void FIgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FIgx2qqxDipoleKernel> documentation
    ("FIgx2qqxDipoleKernel");

}

