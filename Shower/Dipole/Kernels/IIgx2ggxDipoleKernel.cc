// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIgx2ggxDipoleKernel class.
//

#include "IIgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IIgx2ggxDipoleKernel::IIgx2ggxDipoleKernel() 
  : DipoleSplittingKernel() {}

IIgx2ggxDipoleKernel::~IIgx2ggxDipoleKernel() {}

IBPtr IIgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IIgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IIgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool IIgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    a.emitterPDF() == b.emitterPDF() &&
    a.spectatorData() == b.spectatorData() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr IIgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IIgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IIgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IIgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());

  double x = z*(1.-z)/(1.-z+ratio);

  ret *= 3. * ( x/(1.-x) + (1.-x)/x +x*(1.-x) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IIgx2ggxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_mpm = -(1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  //double v_AP_mmp = -v_AP_ppm;
  double v_AP_pmp = -v_AP_mpm;
  
  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_mpm))
    - 2.*abs(rho(0,2))*(v_AP_ppp*v_AP_pmp + v_AP_mmm*v_AP_mpm);

  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*
                                  (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_mpm))/max ) );
  distPhiDep.push_back( make_pair(2, rho(0,2)*(v_AP_mmm*v_AP_mpm + v_AP_pmp*v_AP_ppp)/max ) );
  distPhiDep.push_back( make_pair(-2, rho(2,0)*(v_AP_ppp*v_AP_pmp + v_AP_mpm*v_AP_mmm)/max) );

  return distPhiDep;
}

DecayMEPtr IIgx2ggxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_mpm = -(1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_pmp = -v_AP_mpm;
    
  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1, PDT::Spin1, PDT::Spin1)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm*phase;
  (*kernelPhiDep)(2,2,2) = v_AP_ppp/phase;
  (*kernelPhiDep)(0,0,2) = v_AP_mmp/phase;
  (*kernelPhiDep)(2,2,0) = v_AP_ppm*phase;
  (*kernelPhiDep)(0,2,0) = v_AP_mpm/phase;
  (*kernelPhiDep)(2,0,2) = v_AP_pmp*phase;
  (*kernelPhiDep)(0,2,2) = 0;
  (*kernelPhiDep)(2,0,0) = 0;

  return kernelPhiDep;
}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IIgx2ggxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IIgx2ggxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IIgx2ggxDipoleKernel> IIgx2ggxDipoleKernel::initIIgx2ggxDipoleKernel;
// Definition of the static class description member.

void IIgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<IIgx2ggxDipoleKernel> documentation
    ("IIgx2ggxDipoleKernel");

}

