// -*- C++ -*-
//
// This is a helper header to implement HepMC conversions
//
#include "ThePEG/Vectors/HepMCConverter.h"
#include "HepMC/GenEvent.h"
#include "CLHEP/Vector/LorentzVector.h"

namespace ThePEG {
template<> 
struct HepMCTraits<HepMC::GenEvent> 
  : public HepMCTraitsBase<HepMC::GenEvent,
			   HepMC::GenParticle,
			   HepMC::GenVertex,
			   HepMC::Polarization> 
{
  static ParticleT * newParticle(const Lorentz5Momentum & p,
				 long id, int status) {
    // Note that according to the documentation the momentum is stored in a
    // HepLorentzVector in GeV (event though the CLHEP standard is MeV).
    CLHEP::HepLorentzVector p_GeV(p.x()/GeV, 
				  p.y()/GeV,
				  p.z()/GeV,
				  p.e()/GeV);
    ParticleT * genp = new ParticleT(p_GeV, id, status);
    genp->setGeneratedMass(p.mass()/GeV);
    return genp;
  }

  static void setPosition(VertexT & v, const LorentzPoint & p) {
    // We assume that the position is measured in millimeters.
    CLHEP::HepLorentzVector p_mm(p.x()/mm, 
				 p.y()/mm,
				 p.z()/mm,
				 p.t()/mm);
    v.set_position(p_mm);
  }
};
}
