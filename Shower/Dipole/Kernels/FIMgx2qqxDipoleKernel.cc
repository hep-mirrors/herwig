// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMgx2qqxDipoleKernel class.
//

#include "FIMgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMgx2qqxDipoleKernel::FIMgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

FIMgx2qqxDipoleKernel::~FIMgx2qqxDipoleKernel() {}

IBPtr FIMgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    flavour()->mass() != ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIMgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    emitter(a)->mass() == sk.emitter(b)->mass() &&
    a.spectatorPDF() == b.spectatorPDF();

}


tcPDPtr FIMgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() != ZERO);
  return flavour();
}

tcPDPtr FIMgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6 && flavour()->mass() != ZERO);
  return flavour()->CC();
}

tcPDPtr FIMgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

// TODO
// assure split.scale() is sqrt(sbar)
double FIMgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  // mi=m=mQ, Mi=0, mj=Mj=0
  Energy2 mQ2 = sqr(split.emitterData()->mass());

  double z = split.lastZ();
  double x = 1./ ( 1. +
		   ( sqr(split.lastPt()) + mQ2 ) /
		   ( z*(1.-z) * sqr(split.scale()) ) );

  double muQ2 = x * mQ2/sqr(split.scale());

  double zm = .5 * ( 1. - sqrt( 1. - 4.*muQ2/(1.-x) ) );
  double zp = .5 * ( 1. + sqrt( 1. - 4.*muQ2/(1.-x) ) );

  ret *= .25 * (1.-2.*(zp-z)*(z-zm));

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
FIMgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();
  Energy2 mi2 = sqr(dInfo.emitterData()->mass());
  
  // Altarelli-Parisi spin-indexed kernels:
  double ratio = mi2 / ( mi2 + sqr(pt) );
  double root = sqrt(1.-ratio);
  double v_AP_ppp = sqrt(ratio);
  double v_AP_ppm = z*root;
  double v_AP_pmp = -(1.-z)*root;

  //double v_AP_mmm = v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp)) + 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);
  
  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp) )/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );

  return distPhiDep;
}

DecayMEPtr FIMgx2qqxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {
  
  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();
  Energy2 mi2 = sqr(dInfo.emitterData()->mass());

  // Altarelli-Parisi spin-indexed kernels:
  double ratio = mi2 / ( mi2 + sqr(pt) );
  double root = sqrt(1.-ratio);
  double v_AP_ppp = sqrt(ratio);
  double v_AP_ppm = z*root;
  double v_AP_pmp = -(1.-z)*root;

  double v_AP_mmm = v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  Complex phase = exp(Complex(0.,1.)*dInfo.lastPhi());

  // 0 = -, 2 = +
  (*kernelPhiDep)(0,0,0) = v_AP_mmm;
  (*kernelPhiDep)(2,1,1) = v_AP_ppp;
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


void FIMgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FIMgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIMgx2qqxDipoleKernel> FIMgx2qqxDipoleKernel::initFIMgx2qqxDipoleKernel;
// Definition of the static class description member.

void FIMgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FIMgx2qqxDipoleKernel> documentation
    ("FIMgx2qqxDipoleKernel");

}

