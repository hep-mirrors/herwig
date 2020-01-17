// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFgx2ggxDipoleKernel class.
//

#include "FFgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFgx2ggxDipoleKernel::FFgx2ggxDipoleKernel() 
  : DipoleSplittingKernel(){}

FFgx2ggxDipoleKernel::~FFgx2ggxDipoleKernel() {}

IBPtr FFgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FFgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FFgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

#ifndef NDEBUG
bool FFgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
#else
bool FFgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& ,
#endif
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g;
       

}

tcPDPtr FFgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FFgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FFgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double y = sqr(split.lastPt() / split.scale()) / (z*(1.-z));

  double S1=1./(1.-z*(1.-y));
  double S2=1./(1.-(1.-z)*(1.-y));
  double NS=(-2 + z*(1.-z));
  
  if( theAsymmetryOption == 0 ){
    ret *= 3.*( S1 + 0.5 * NS);
  }else if ( theAsymmetryOption == 1 ){
    ret *= 3.*z*( S1 +S2 + NS );
  }else{
    ret *= 3.*0.5*( S1 + S2 + NS );
  }
  
  return ret > 0. ? ret : 0.;
}

vector< pair<int, Complex> >
FFgx2ggxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const {

  double z = dInfo.lastZ();

  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_pmp = (1.-z)*sqrt( (1.-z)/z );

  //double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Initialise variables for the distributions
  vector< pair<int, Complex> > distPhiDep;
  double max = (sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp)) - 2.*abs(rho(0,2))*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp);

  distPhiDep.push_back( make_pair(0, (rho(0,0)+rho(2,2))*(sqr(v_AP_ppp) + sqr(v_AP_ppm) + sqr(v_AP_pmp))/max ) );
  distPhiDep.push_back( make_pair(-2, rho(0,2)*(v_AP_mpm*v_AP_ppm + v_AP_mmp*v_AP_pmp)/max ) );
  distPhiDep.push_back( make_pair(2, rho(2,0)*(v_AP_ppm*v_AP_mpm + v_AP_pmp*v_AP_mmp)/max) );

  return distPhiDep;
}

DecayMEPtr FFgx2ggxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  double z = dInfo.lastZ();
  
  // Altarelli-Parisi spin-indexed kernels:
  double v_AP_ppp = -sqrt( 1./(z*(1.-z)) );
  double v_AP_ppm = z*sqrt( z / (1.-z) );
  double v_AP_pmp = (1.-z)*sqrt( (1.-z)/z );

  double v_AP_mmm = -v_AP_ppp;
  double v_AP_mmp = -v_AP_ppm;
  double v_AP_mpm = -v_AP_pmp;

  // Construct the (phi-dependent) spin-unaveraged splitting kernel
  DecayMEPtr kernelPhiDep
    (new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
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


void FFgx2ggxDipoleKernel::persistentOutput(PersistentOStream & os) const {
  os << theAsymmetryOption;
}

void FFgx2ggxDipoleKernel::persistentInput(PersistentIStream & is, int) {
  is >> theAsymmetryOption;
}

ClassDescription<FFgx2ggxDipoleKernel> FFgx2ggxDipoleKernel::initFFgx2ggxDipoleKernel;
// Definition of the static class description member.

void FFgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<FFgx2ggxDipoleKernel> documentation
    ("FFgx2ggxDipoleKernel");

  static Parameter<FFgx2ggxDipoleKernel,int> interfacetheAsymmetryOption
    ("AsymmetryOption",
     "The asymmetry option for final state gluon spliitings.",
     &FFgx2ggxDipoleKernel::theAsymmetryOption, 1, 0, 0,
     false, false, Interface::lowerlim);
}
