// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AnomalousVVVVertex class.
//

#include "AnomalousVVVVertex.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace RadiativeZPrime;

AnomalousVVVVertex::AnomalousVVVVertex() {
  setNpoint(3);
  setSpin(3,3,3);
  setName(VVV);
}

AbstractNoPIOClassDescription<AnomalousVVVVertex> 
AnomalousVVVVertex::initAnomalousVVVVertex;
// Definition of the static class description member.

void AnomalousVVVVertex::Init() {

  static ClassDocumentation<AnomalousVVVVertex> documentation
    ("There is no documentation for the AnomalousVVVVertex class");

}

// evaluate the vertex
Complex AnomalousVVVVertex::evaluate(Energy2 q2, const VectorWaveFunction & vec1,
			    const VectorWaveFunction & vec2,
			    const VectorWaveFunction & vec3) {
  // calculate the coupling
  setCoupling(q2,vec1.getParticle(),vec2.getParticle(),vec3.getParticle());
  LorentzPolarizationVector eps = epsilon(vec1.wave(),vec2.wave(),vec3.wave());
  if(vec1.getParticle()->id()==ParticleID::gamma) 
    return getNorm()*Complex(0.,1.)*(eps*vec1.getMomentum())*UnitRemoval::InvE;
  else if(vec2.getParticle()->id()==ParticleID::gamma) 
    return getNorm()*Complex(0.,1.)*(eps*vec2.getMomentum())*UnitRemoval::InvE;
  else
    return getNorm()*Complex(0.,1.)*(eps*vec3.getMomentum())*UnitRemoval::InvE;
}

// off-shell vector
VectorWaveFunction AnomalousVVVVertex::evaluate(Energy2, int, tcPDPtr,
						const VectorWaveFunction & ,
						const VectorWaveFunction &,
						Energy, Energy) {
  throw Exception() << "AnomalousVVVVertex::evaluate() only implemented "
		    << "for the evaluation of the vertex, not for the "
		    << "evaluation of the off-shell vector wavefunction"
		    << Exception::runerror;
}
