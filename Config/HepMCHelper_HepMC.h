// -*- C++ -*-
//
// This is a helper header to implement HepMC conversions
//
#include "ThePEG/Vectors/HepMCConverter.h"
#include "HepMC/GenEvent.h"

namespace ThePEG {
template<> 
struct HepMCTraits<HepMC::GenEvent> 
  : public HepMCTraitsBase<HepMC::GenEvent,
			   HepMC::GenParticle,
			   HepMC::GenVertex,
			   HepMC::Polarization> 
{};
}
