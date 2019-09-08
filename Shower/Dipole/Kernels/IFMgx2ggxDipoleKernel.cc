// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMgx2ggxDipoleKernel class.
//

#include "IFMgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMgx2ggxDipoleKernel::IFMgx2ggxDipoleKernel() 
  : DipoleSplittingKernel() {}

IFMgx2ggxDipoleKernel::~IFMgx2ggxDipoleKernel() {}

IBPtr IFMgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFMgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFMgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() != ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFMgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    a.emitterPDF() == b.emitterPDF() &&
    sk.spectator(b)->mass() == spectator(a)->mass();

}


tcPDPtr IFMgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFMgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFMgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFMgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  Energy pt = split.lastPt();
  double ratio = sqr(pt/split.scale());
  double muk2 = sqr(split.spectatorMass()/split.scale());

	// Calculate x and u
    double rho = 1. - 4.*ratio*(1.-muk2)*z*(1.-z)/sqr(1.-z+ratio);
    double x = 0.5*((1.-z+ratio)/(ratio*(1.-muk2))) * (1. - sqrt(rho));
    double u = x*ratio / (1.-z);

// NOTE - The definition of muk used in the kinematics differs from that in CS

    double muk2CS = x*muk2;
    ret *= 3. * ( 1./(1.-x+u) - 1. + x*(1.-x) + (1.-x)/x - muk2CS*u/(x*(1.-u)) );

  return ret > 0. ? ret : 0.;

}

vector< pair<int, Complex> >
IFMgx2ggxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

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
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_mpm)) - 2.*abs(rho(0,2))*(v_AP_ppp*v_AP_pmp + v_AP_mmm*v_AP_mpm);

  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_mpm))/max ) );
  distPhiDep.push_back( make_pair(2, rho(0,2)*(v_AP_mmm*v_AP_mpm + v_AP_pmp*v_AP_ppp)/max ) );
  distPhiDep.push_back( make_pair(-2, rho(2,0)*(v_AP_ppp*v_AP_pmp + v_AP_mpm*v_AP_mmm)/max) );
  
  return distPhiDep;
}

DecayMEPtr IFMgx2ggxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

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


void IFMgx2ggxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFMgx2ggxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFMgx2ggxDipoleKernel> IFMgx2ggxDipoleKernel::initIFMgx2ggxDipoleKernel;
// Definition of the static class description member.

void IFMgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<IFMgx2ggxDipoleKernel> documentation
    ("IFMgx2ggxDipoleKernelv");

}

