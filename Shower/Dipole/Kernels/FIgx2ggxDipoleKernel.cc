// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIgx2ggxDipoleKernel class.
//

#include "FIgx2ggxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIgx2ggxDipoleKernel::FIgx2ggxDipoleKernel() 
  : DipoleSplittingKernel() {}

FIgx2ggxDipoleKernel::~FIgx2ggxDipoleKernel() {}

IBPtr FIgx2ggxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr FIgx2ggxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool FIgx2ggxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
  useThisKernel() &&
    ind.emitterData()->id() == ParticleID::g &&
    ind.spectatorData()->mass() == ZERO &&
    !ind.initialStateEmitter() && ind.initialStateSpectator();
}

bool FIgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    sk.emitter(b)->id() == ParticleID::g &&
    sk.emission(b)->id() == ParticleID::g &&
    a.spectatorPDF() == b.spectatorPDF();

}

tcPDPtr FIgx2ggxDipoleKernel::emitter(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIgx2ggxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr FIgx2ggxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double FIgx2ggxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double x = 1. / ( 1. + sqr(split.lastPt()/split.scale()) / (z*(1.-z)) );

  double S1=1./(1.-z+(1.-x));
  double S2=1./(1.-(1.-z)+(1.-x));
  double NS=(-2 + z*(1.-z)+(1.-x)*(1.+x*z*(1.-z)));
  
  if( theAsymmetryOption == 0 ){
    ret *= 3.*( S1 + 0.5 * NS);
  }else if ( theAsymmetryOption == 1 ){
    ret *= 3.*z*( S1 +S2 + NS );
  }else{
    ret *= 3.*0.5*( S1 + S2 + NS );
  }
  
  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIgx2ggxDipoleKernel::persistentOutput(PersistentOStream & os) const {

  os << theAsymmetryOption;
}

void FIgx2ggxDipoleKernel::persistentInput(PersistentIStream & is, int) {
  is >> theAsymmetryOption;
}

ClassDescription<FIgx2ggxDipoleKernel> FIgx2ggxDipoleKernel::initFIgx2ggxDipoleKernel;
// Definition of the static class description member.

void FIgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<FIgx2ggxDipoleKernel> documentation
    ("FIgx2ggxDipoleKernel");

  static Parameter<FIgx2ggxDipoleKernel,int> interfacetheAsymmetryOption
  ("AsymmetryOption",
   "The asymmetry option for final state gluon spliitings.",
   &FIgx2ggxDipoleKernel::theAsymmetryOption, 1, 0, 0,
   false, false, Interface::lowerlim);
}

