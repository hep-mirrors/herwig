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
				 long id, int status, Energy unit) {
    CLHEP::HepLorentzVector p_unit(p.x()/unit, 
				   p.y()/unit,
				   p.z()/unit,
				   p.e()/unit);
    ParticleT * genp = new ParticleT(p_unit, id, status);
    genp->setGeneratedMass(p.mass()/unit);
    return genp;
  }

  /**
   *  Set the position
   */
  static void setPosition(VertexT & v, const LorentzPoint & p, 
			  Length unit) {
    // We assume that the position is measured in millimeters.
    CLHEP::HepLorentzVector p_unit(p.x()/unit, 
				   p.y()/unit,
				   p.z()/unit,
				   p.t()/unit);
    v.set_position(p_unit);
  }
};
}
