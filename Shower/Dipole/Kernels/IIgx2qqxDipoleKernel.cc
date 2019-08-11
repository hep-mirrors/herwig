// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIgx2qqxDipoleKernel class.
//

#include "IIgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IIgx2qqxDipoleKernel::IIgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

IIgx2qqxDipoleKernel::~IIgx2qqxDipoleKernel() {}

IBPtr IIgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IIgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IIgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() == ZERO &&
    ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool IIgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    flavour() == sk.flavour() &&
    abs(flavour()->id()) < 6 &&
    flavour()->mass() == ZERO &&
    a.emitterPDF() == b.emitterPDF() &&
    a.spectatorData() == b.spectatorData() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr IIgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IIgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() == ZERO);
  return flavour();
}

tcPDPtr IIgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IIgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());

  double x = z*(1.-z)/(1.-z+ratio);

  ret *= 0.5 * (!strictLargeN() ? 4./3. : 3./2.) * ( 1./x +sqr(1.-x)/x );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IIgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {
  
  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt(1./z);
  double v_AP_mpm = (1.-z)/sqrt(z);

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_pmp = -v_AP_mpm;
    
  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;  
  double max = sqr(v_AP_ppp) + sqr(v_AP_mpm)
    - 2.*abs(rho(0,2))*(v_AP_pmp*v_AP_ppp + v_AP_mmm*v_AP_mpm);

  distPhiDep.push_back(make_pair( 0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_mpm))/max ) );
  distPhiDep.push_back(make_pair( 2, rho(0,2)*(v_AP_pmp*v_AP_ppp + v_AP_mmm*v_AP_mpm )/max ) );
  distPhiDep.push_back(make_pair( -2, rho(2,0)*(v_AP_ppp*v_AP_pmp + v_AP_mpm*v_AP_mmm )/max ) );

  return distPhiDep;
}

DecayMEPtr IIgx2qqxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = sqrt(1./z);
  double v_AP_mpm = (1.-z)/sqrt(z);
  
  double v_AP_mmm = -v_AP_ppp;
  double v_AP_pmp = -v_AP_mpm;
  
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half, PDT::Spin1, PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -1 or -1/2, 1=+1/2, 2 = +1
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(1,2,1) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,1) = 0.;
  (*kernelPhiDep)(1,2,0) = 0.;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(1,0,1) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,1) = 0.;
  (*kernelPhiDep)(1,0,0) = 0.;            
  
  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IIgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IIgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IIgx2qqxDipoleKernel> IIgx2qqxDipoleKernel::initIIgx2qqxDipoleKernel;
// Definition of the static class description member.

void IIgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<IIgx2qqxDipoleKernel> documentation
    ("IIgx2qqxDipoleKernel");

}

