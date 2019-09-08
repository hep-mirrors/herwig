// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMgx2qqxDipoleKernel class.
//

#include "FFMgx2qqxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFMgx2qqxDipoleKernel::FFMgx2qqxDipoleKernel() 
  : DipoleSplittingKernel() {}

FFMgx2qqxDipoleKernel::~FFMgx2qqxDipoleKernel() {}

IBPtr FFMgx2qqxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFMgx2qqxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFMgx2qqxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    !( ind.spectatorData()->mass() == ZERO &&
       flavour()->mass() == ZERO ) &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator() &&
    !ind.incomingDecayEmitter() && !ind.incomingDecaySpectator();
}

bool FFMgx2qqxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {
  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() + sk.emission(b)->id() == 0 &&
    abs(sk.emitter(b)->id()) < 6 &&
    emitter(a)->id() == sk.emitter(b)->id() &&
    abs(sk.spectator(b)->mass()) == abs(spectator(a)->mass());
}

tcPDPtr FFMgx2qqxDipoleKernel::emitter(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6);
  return flavour();
}

tcPDPtr FFMgx2qqxDipoleKernel::emission(const DipoleIndex&) const {
  assert(flavour());
  assert(abs(flavour()->id()) < 6);
  return flavour()->CC();
}

tcPDPtr FFMgx2qqxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFMgx2qqxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);

  // Sudakov parameterisation variables,
  // needed to calculate y.
  double z = split.lastZ();
  Energy pt = split.lastPt();

  // Construct mass squared variables
  Energy2 Qijk = sqr(split.scale());
  Energy2 mi2 = sqr(split.emitterData()->mass());
  Energy2 mj2 = mi2;
  Energy2 mk2 = sqr(split.spectatorMass());
  Energy2 sbar = Qijk - mi2 - mj2 - mk2;

  // Calculate y
  double y = (sqr(pt) + sqr(1.-z)*mi2 + sqr(z)*mj2) / sbar / z / (1.-z);
  
  // zi, used in dipole splitting kernel
  double zi = split.lastSplittingParameters()[0];
  
  double vijk = sqrt( sqr(2.*mk2 + sbar*(1.-y)) - 4.*mk2*Qijk ) / sbar / (1.-y);
  double viji = sqrt( sqr(sbar*y) - 4.*sqr(mi2) ) / (sbar*y + 2.*mi2);

  double zip = 0.5*(1.+viji*vijk);
  double zim = 0.5*(1.-viji*vijk);

  // how to choose kappa??
  double kappa = 0.;

  ret *= 0.25 / vijk *
    ( 1. - 2.*( zi*(1.-zi) - (1.-kappa)*zip*zim 
                - kappa*mi2 / ( 2.*mi2 + (Qijk - 2.*mi2 - mk2)*y) ) );
    
  return ret > 0. ? ret : 0.;
  
}

vector< pair<int, Complex> >
FFMgx2qqxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  // Need variables for the AP kernels
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

DecayMEPtr FFMgx2qqxDipoleKernel::matrixElement( const DipoleSplittingInfo& dInfo ) const {
  
  // Need variables for the AP kernels
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


void FFMgx2qqxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void FFMgx2qqxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFMgx2qqxDipoleKernel> 
FFMgx2qqxDipoleKernel::initFFMgx2qqxDipoleKernel;
// Definition of the static class description member.

void FFMgx2qqxDipoleKernel::Init() {

  static ClassDocumentation<FFMgx2qqxDipoleKernel> documentation
    ("FFMgx2qqxDipoleKernel");

}

