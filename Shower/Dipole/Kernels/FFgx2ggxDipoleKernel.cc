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
  : DipoleSplittingKernel(),theSymmetryFactor(0.5){}

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

  ret *=theSymmetryFactor*3.*(1./(1.-z*(1.-y))+1./(1.-(1.-z)*(1.-y))-2.+z*(1.-z));
  

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFgx2ggxDipoleKernel::persistentOutput(PersistentOStream & os) const {

  os<<theSymmetryFactor;
}

void FFgx2ggxDipoleKernel::persistentInput(PersistentIStream & is, int) {
  is>>theSymmetryFactor;
}

ClassDescription<FFgx2ggxDipoleKernel> FFgx2ggxDipoleKernel::initFFgx2ggxDipoleKernel;
// Definition of the static class description member.

void FFgx2ggxDipoleKernel::Init() {

  static ClassDocumentation<FFgx2ggxDipoleKernel> documentation
    ("FFgx2ggxDipoleKernel");

  static Parameter<FFgx2ggxDipoleKernel,double> interfaceSymmetryFactor
    ("SymmetryFactor",
     "The symmetry factor for final state gluon spliitings.",
     &FFgx2ggxDipoleKernel::theSymmetryFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);
}

