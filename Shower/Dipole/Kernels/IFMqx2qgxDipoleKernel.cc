// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMqx2qgxDipoleKernel class.
//

#include "IFMqx2qgxDipoleKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMqx2qgxDipoleKernel::IFMqx2qgxDipoleKernel() 
  : DipoleSplittingKernel() {}

IFMqx2qgxDipoleKernel::~IFMqx2qgxDipoleKernel() {}

IBPtr IFMqx2qgxDipoleKernel::clone() const {
  return new_ptr(*this);
}

IBPtr IFMqx2qgxDipoleKernel::fullclone() const {
  return new_ptr(*this);
}

bool IFMqx2qgxDipoleKernel::canHandle(const DipoleIndex& ind) const {
  return
    abs(ind.emitterData()->id()) < 6  &&
    ind.emitterData()->mass() == ZERO &&
    ind.spectatorData()->mass() != ZERO &&
    ind.initialStateEmitter() && !ind.initialStateSpectator();
}

bool IFMqx2qgxDipoleKernel::canHandleEquivalent(const DipoleIndex& a,
					       const DipoleSplittingKernel& sk,
					       const DipoleIndex& b) const {

  assert(canHandle(a));

  if ( !canHandle(b) )
    return false;

  return
    emitter(a) == sk.emitter(b) &&
    emission(a) == sk.emission(b) &&
    spectator(a)->mass() == sk.spectator(b)->mass() &&
    a.emitterPDF() == b.emitterPDF();

}


tcPDPtr IFMqx2qgxDipoleKernel::emitter(const DipoleIndex& ind) const {
  return ind.emitterData();
}

tcPDPtr IFMqx2qgxDipoleKernel::emission(const DipoleIndex&) const {
  return getParticleData(ParticleID::g);
}

tcPDPtr IFMqx2qgxDipoleKernel::spectator(const DipoleIndex& ind) const {
  return ind.spectatorData();
}

double IFMqx2qgxDipoleKernel::evaluate(const DipoleSplittingInfo& split) const {

  double ret = alphaPDF(split);

  double z = split.lastZ();
  double ratio = sqr(split.lastPt()/split.scale());
  double alpha = 1. - 2.*sqr(split.spectatorData()->mass()/split.scale());
  double root  = sqr(1.-z+alpha*ratio) - 4.*ratio*(1.-z);
  if(root < 0. && root > -1e-10)
    root = 0.;
  else if(root < 0.)
    return 0.;
  root = sqrt(root);

  double x = alpha == 1. ? ( z*(1.-z) - ratio ) / ( 1. - z - ratio ) :
    ( sqr(alpha)*ratio + 2.*z - alpha*(1.+z) + alpha*root ) /
    (2.*(1.-alpha));
  double u = ( 1.-z + alpha*ratio - root ) /
    (2.*(1.-z));

  ret *= (!strictLargeN() ? 4./3. : 3./2.) * ( 2./(1.-x+u) - (1.+x) + u*(1.+3.*x*(1.-u) ) );

  return ret > 0. ? ret : 0.;

}

// If needed, insert default implementations of  function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMqx2qgxDipoleKernel::persistentOutput(PersistentOStream & ) const {
}

void IFMqx2qgxDipoleKernel::persistentInput(PersistentIStream & , int) {
}

ClassDescription<IFMqx2qgxDipoleKernel> IFMqx2qgxDipoleKernel::initIFMqx2qgxDipoleKernel;
// Definition of the static class description member.

void IFMqx2qgxDipoleKernel::Init() {

  static ClassDocumentation<IFMqx2qgxDipoleKernel> documentation
    ("IFMqx2qgxDipoleKernel");

}

