// -*- C++ -*-
//
// HepMCHelper_HepMC.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
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

namespace ThePEG {
/**
 * Struct for HepMC conversion
 */
template<> 
struct HepMCTraits<HepMC::GenEvent> 
  : public HepMCTraitsBase<HepMC::GenEvent,
			   HepMC::GenParticle,
			   HepMC::GenVertex,
			   HepMC::Polarization,
			   HepMC::PdfInfo>
{};
}
