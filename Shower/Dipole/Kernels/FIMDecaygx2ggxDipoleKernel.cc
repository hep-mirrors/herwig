// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMDecaygx2ggxDipoleKernel class.
//

#include "FIMDecaygx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMDecaygx2ggxDipoleKernel::FIMDecaygx2ggxDipoleKernel() 
  : DipoleSplittingKernel(){}

FIMDecaygx2ggxDipoleKernel::~FIMDecaygx2ggxDipoleKernel() {}

IBPtr FIMDecaygx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIMDecaygx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIMDecaygx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.incomingDecaySpectator() && !ind.incomingDecayEmitter() &&
    ind.emitterData()->id() == ParticleID::g &&
    !(ind.spectatorData()->mass() == ZERO) &&
    // Initial state here refers to the entire event
    !ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool FIMDecaygx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emission(b)->id() == ParticleID::g &&
    sk.emitter(b)->id() == ParticleID::g &&
    abs(sk.spectator(b)->mass()) == abs(spectator(a)->mass());

}


tcPDPtr FIMDecaygx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIMDecaygx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIMDecaygx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIMDecaygx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {
  
  double ret = alphaPDF(split);

  // Sudakov parameterisation variables,
  // needed to calculate y.
  double z = split.lastZ();
  Energy pt = split.lastPt();

  // Construct mass squared variables
  // Note for q->qg can use the emitterMass
  // (i.e. mass of emitter before splitting = mass of emitter after)
  Energy2 Qijk = sqr(split.scale());
  Energy2 mk2 = sqr(split.recoilMass());
  Energy2 sbar = Qijk - mk2;

  // Note this should be the same as Qijk
  Energy2 ma2 = sqr(split.spectatorMass());


  // Calculate y
  double y = sqr(pt) / sbar / z / (1.-z);

  if( sqr(2.*mk2+sbar*(1.-y)) - 4.*mk2*Qijk < ZERO ){
    generator()->logWarning( Exception()
    << "error in FIMDecayqx2qgxDipoleKernel::evaluate -- " <<
    "mk2 " << mk2/GeV2 << "  y " << y << Exception::warning );
    return 0.0;
  }

  // zi, used in dipole splitting kernel
  double zi = split.lastSplittingParameters()[0];
  
  double vijk = sqrt( sqr(2.*mk2 + sbar*(1.-y)) - 4.*mk2*Qijk ) / sbar / (1.-y);
  double vtilde = 1.;
  double viji = 1.;

  double zip = 0.5*(1.+viji*vijk);
  double zim = 0.5*(1.-viji*vijk);

  // how to choose kappa?
  double kappa = 0.;

  double S1 = 0.5*3.*(2.*y + 1.)/((1.+y)-zi*(1.-y)) +
    (!strictLargeN() ? 4./3. : 3./2.)*
    y/(1.-zi*(1.-y)) * ( 2.*(2.*y + 1.)/((1.+y)-zi*(1.-y))
			- (vtilde/vijk)*(2. + 2.*ma2/((1.-zi*(1.-y))*sbar)) );
  double S2 = 0.5*3.*(2.*y + 1.)/((1.+y)-(1.-zi)*(1.-y)) +
    (!strictLargeN() ? 4./3. : 3./2.)*
    y/(1.-(1.-zi)*(1.-y)) * ( 2.*(2.*y + 1.)/((1.+y)-(1.-zi)*(1.-y))
			     - (vtilde/vijk)*(2. + 2.*ma2/((1.-(1.-zi)*(1.-y))*sbar)) );  
  double NS = 0.5*3.*(zi*(1.-zi)-(1.-kappa)*zip*zim - 2.)/vijk;

  if( theAsymmetryOption == 0 ){
    ret *= 2.*S1 + NS;
  }else if ( theAsymmetryOption == 1 ){
    ret *= 2.*zi*( S1 + S2 + NS );
  }else{
    ret *= S1 + S2 + NS;
  }
  
  return ret > 0. ? ret : 0.;
  
}

vector< pair<int, Complex> >
FIMDecaygx2ggxDipoleKernel::generatePhi(const DipoleSplittingInfo& dInfo,
                                        const RhoDMatrix& rho) const {

  // Need variables for the AP kernels
  double z = dInfo.lastSplittingParameters()[0];

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


DecayMEPtr FIMDecaygx2ggxDipoleKernel::matrixElement(const DipoleSplittingInfo& dInfo) const {

  // Need variables for the AP kernels
  double z = dInfo.lastSplittingParameters()[0];  
  
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


void FIMDecaygx2ggxDipoleKernel::persistentOutput(PersistentOStream & os) const {
  os<<theAsymmetryOption;
}

void FIMDecaygx2ggxDipoleKernel::persistentInput(PersistentIStream & is, int) {
  is>>theAsymmetryOption;
}

ClassDescription<FIMDecaygx2ggxDipoleKernel> FIMDecaygx2ggxDipoleKernel::initFIMDecaygx2ggxDipoleKernel;
// Definition of the static class description member.

void FIMDecaygx2ggxDipoleKernel::Init() {

  static ClassDocumentation<FIMDecaygx2ggxDipoleKernel> documentation
    ("FIMDecaygx2ggxDipoleKernel");

  static Parameter<FIMDecaygx2ggxDipoleKernel,int> interfacetheAsymmetryOption
    ("AsymmetryOption",
     "The asymmetry option for final state gluon spliitings.",
     &FIMDecaygx2ggxDipoleKernel::theAsymmetryOption, 0, 0, 0,
     false, false, Interface::lowerlim);

}
