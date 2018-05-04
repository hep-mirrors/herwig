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

bool FFgx2ggxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
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
