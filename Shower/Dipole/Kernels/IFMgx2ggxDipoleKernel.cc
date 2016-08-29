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
  double ratio = sqr(split.lastPt()/split.scale());
  double muj2 = sqr(split.spectatorData()->mass()/split.scale());
  double alpha = 1. - 2.*muj2;
  double root  = sqr(1.-z+alpha*ratio) - 4.*ratio*(1.-z);
  if(root < 0. && root > -1e-10)
    root = 0.;
  else if(root < 0.)
    return 0.;
  root = sqrt(root);

  double x = ( sqr(alpha)*ratio + 2.*z - alpha*(1.+z) + alpha*root ) /
    (2.*(1.-alpha));
  double u = ( 1.-z + alpha*ratio - root ) /
    (2.*(1.-z));

  // careful: CSmassless u is CSmassive z_i
  ret *= 3. * ( 1./(1.-x+u) - 1. + x*(1.-x) + (1.-x)/x - muj2*u/(x*(1.-u)) );

  return ret > 0. ? ret : 0.;

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

