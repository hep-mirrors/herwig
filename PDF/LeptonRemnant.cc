// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptonRemnant class.
//

#include "LeptonRemnant.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/UseRandom.h"

using namespace Herwig;

void LeptonRemnant::doinit() throw(InitException) {
  photon = getParticleData(ParticleID::gamma);
  RemnantHandler::doinit();
}

void LeptonRemnant::persistentOutput(PersistentOStream & os) const {
  os << photon << minX;
}

void LeptonRemnant::persistentInput(PersistentIStream & is, int) {
  is >> photon >> minX;
}

ClassDescription<LeptonRemnant> LeptonRemnant::initLeptonRemnant;
// Definition of the static class description member.

void LeptonRemnant::Init() {

  static ClassDocumentation<LeptonRemnant> documentation
    ("The LeptonRemnant class is design to generate the remnant in leptonic collisions");

  static Parameter<LeptonRemnant,double> interfaceMinX
    ("MinX",
     "The minimum energy fraction allowed for a photon remnant. "
     "If less than this no remnant will be emitted.",
     &LeptonRemnant::minX, 1.0e-10, 0.0, 1.0,
     true, false, true);
}

bool LeptonRemnant::
canHandle(tcPDPtr particle, const cPDVector & partons) const {
  // check that we have an incoming lepton beam
  int id = abs(particle->id());
  if( id != ParticleID::eminus && id != ParticleID::muminus) return false;
  for ( cPDVector::const_iterator it = partons.begin(); it != partons.end(); ++it ) {
    if ( (**it).id() != particle->id()    &&
	 (**it).id() != ParticleID::gamma ) return false;
  }
  return true;
}

Lorentz5Momentum LeptonRemnant::
generate(PartonBinInstance & pb, const double *,
	 Energy2 scale, const LorentzMomentum & parent) const {
  // photon into hard process and lepton remnant
  if ( pb.particleData() != pb.partonData() && 
       pb.partonData()->id() == ParticleID::gamma) {
    Energy  ppl = pb.xi()*(abs(parent.z())+parent.t());
    Energy2 qt2 = pb.eps()*scale-sqr(pb.xi()*parent.m());
    Energy  pmi = (qt2-scale)/ppl;
    Lorentz5Momentum pgam;
    pgam.setMass(-sqrt(scale));
    pgam.setT(0.5*(ppl+pmi));
    pgam.setZ(0.5*(ppl-pmi));
    if(parent.z()<0.*GeV) pgam.z() *=-1.;
    double phi = rnd(2.0*Constants::pi);
    pgam.setX(sqrt(qt2)*cos(phi));
    pgam.setY(sqrt(qt2)*sin(phi));
    pgam.rotateY(parent.theta());
    pgam.rotateZ(parent.phi());
    Lorentz5Momentum prem=parent-pgam;
    PPtr rem = pb.particleData()->produceParticle(prem, pb.particleData()->mass());
    pb.remnants(PVector(1, rem));
    return pgam;
  }
  else if( pb.particleData() == pb.partonData() ) {
    if ( pb.eps() < minX ) {
      pb.remnants(PVector());
      return parent;
    }
    LorentzMomentum p(0.0*GeV, 0.0*GeV, parent.rho(), parent.e());
    TransverseMomentum qt;
    Energy2 qt2 = 0.0*GeV2;
    if ( scale >= 0.0*GeV2 ) {
      qt2 = pb.eps()*(pb.xi()*parent.m2() + scale);
      double phi = rnd(2.0*Constants::pi);
      qt = TransverseMomentum(sqrt(qt2)*cos(phi), sqrt(qt2)*sin(phi));
    }
    Energy pl = p.plus()*pb.eps();
    LorentzMomentum prem = lightCone(pl, qt2/pl, qt);
    prem.rotateY(parent.theta());
    prem.rotateZ(parent.phi());
    PPtr rem = photon->produceParticle(prem, 0.0*GeV);
    pb.remnants(PVector(1, rem));
    return parent - rem->momentum();
  }
  else {
    throw RemnantHandlerException
      (pb.particleData()->name(), pb.partonData()->name(), name(),
       "The remnant handler can only extract leptons of the"
       " same type and photons from leptons.");
  }
}

Lorentz5Momentum LeptonRemnant::
generate(PartonBinInstance & pb, const double *, Energy2 scale, Energy2,
	 const LorentzMomentum & parent) const {
  // photon into hard process and lepton remnant
  if ( pb.particleData() != pb.partonData() && 
       pb.partonData()->id() == ParticleID::gamma) {
    Energy  ppl = pb.xi()*(abs(parent.z())+parent.t());
    Energy2 qt2 = pb.eps()*scale-sqr(pb.xi()*parent.m());
    Energy  pmi = (qt2-scale)/ppl;
    Lorentz5Momentum pgam;
    pgam.setMass(-sqrt(scale));
    pgam.setT(0.5*(ppl+pmi));
    pgam.setZ(0.5*(ppl-pmi));
    if(parent.z()<0.*GeV) pgam.z() *=-1.;
    double phi = rnd(2.0*Constants::pi);
    pgam.setX(sqrt(qt2)*cos(phi));
    pgam.setY(sqrt(qt2)*sin(phi));
    pgam.rotateY(parent.theta());
    pgam.rotateZ(parent.phi());
    Lorentz5Momentum prem=parent-pgam;
    PPtr rem = pb.particleData()->produceParticle(prem, pb.particleData()->mass());
    pb.remnants(PVector(1, rem));
    return pgam;
  }
  else if ( pb.particleData() == pb.partonData() ) {
    if ( pb.eps() < minX ) {
      pb.remnants(PVector());
      return parent;
    }
    LorentzMomentum p(0.0*GeV, 0.0*GeV, parent.rho(), parent.e());
    TransverseMomentum qt;
    Energy2 qt2 = 0.0*GeV2;
    if ( scale >= 0.0*GeV2 ) {
      qt2 = pb.eps()*(pb.xi()*parent.m2() + scale);
      double phi = rnd(2.0*Constants::pi);
      qt = TransverseMomentum(sqrt(qt2)*cos(phi), sqrt(qt2)*sin(phi));
    }
    Energy pl = p.plus()*pb.eps();
    LorentzMomentum prem = lightCone(pl, qt2/pl, qt);
    prem.rotateY(parent.theta());
    prem.rotateZ(parent.phi());
    PPtr rem = photon->produceParticle(prem, 0.0*GeV);
    pb.remnants(PVector(1, rem));
    return parent - rem->momentum();
  }
  else {
    throw RemnantHandlerException
      (pb.particleData()->name(), pb.partonData()->name(), name(),
       "The remnant handler can only extract leptons of the"
       " same type and photons from leptons.");
  }
}
