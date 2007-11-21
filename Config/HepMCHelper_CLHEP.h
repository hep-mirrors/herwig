// -*- C++ -*-
//
// HepMCHelper_CLHEP.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
  /**
   *  Create a new particle
   * @param p The momentum
   * @param id The id
   * @param status The status
   */
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

  /**
   *  Set the position
   */
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
