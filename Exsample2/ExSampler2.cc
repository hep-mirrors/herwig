// -*- C++ -*-
//
// ExSampler.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ExSampler class.
//

#include "ExSampler2.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "GeneralSampler.h"

using namespace Herwig;
using namespace exsample;

ExSampler::ExSampler() 
  : presampling_points_(1000),
    freeze_grid_(0),
    efficiency_threshold_(.95),
    gain_threshold_(.1) {}

ExSampler::~ExSampler() {}

IBPtr ExSampler::clone() const {
  return new_ptr(*this);
}

IBPtr ExSampler::fullclone() const {
  return new_ptr(*this);
}

void ExSampler::generate(bool) {

  while (true) {
    try {
      generator_.generate(*this);
      lastPoint() = generator_.last_point();
      break;

    } catch(selection_maxtry&) {
      throw GeneralSampler::MaxTryException()
	<< "The maximum number of attempts to select a cell was exceeded\n"
	<< "for process " << process()
	<< Exception::eventerror;
    } catch(hit_and_miss_maxtry&) {
      throw GeneralSampler::MaxTryException()
	<< "The maximum number of attempts to select an event was exceeded\n"
	<< "for process " << process()
	<< Exception::eventerror;
    } catch(generator_update&) {
      throw UpdateCrossSections();
    } catch(...) {
      throw;
    }
  }

}

void ExSampler::initialize(bool progress) {

  if ( progress ) {
    cout << "initializing sampler for "
	 << process() << "\n";
  }

  generator_.function(this);

  generator_.sampling_parameters().presampling_points = presampling_points_;
  generator_.sampling_parameters().freeze_grid = freeze_grid_;
  generator_.sampling_parameters().maxtry = eventHandler()->maxLoop();
  generator_.sampling_parameters().efficiency_threshold = efficiency_threshold_;
  generator_.sampling_parameters().gain_threshold = gain_threshold_;

  generator_.initialize(*this);

  if ( progress ) {
    cout << "estimated cross section is ( "
	 << averageWeight() << " +/- " << sqrt(averageWeightVariance())
	 << " ) nb\n";
  }

}

void ExSampler::finalize(bool) {
  generator_.finalize();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void ExSampler::persistentOutput(PersistentOStream & os) const {
  os << presampling_points_ << freeze_grid_ 
     << efficiency_threshold_ << gain_threshold_;
  generator_.put(os);
}

void ExSampler::persistentInput(PersistentIStream & is, int) {
  is >> presampling_points_ >> freeze_grid_ 
     >> efficiency_threshold_ >> gain_threshold_;
  generator_.get(is);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ExSampler,Herwig::BinSampler>
  describeHerwigExSampler("Herwig::ExSampler", "HwExsample2.so");

void ExSampler::Init() {

  static ClassDocumentation<ExSampler> documentation
    ("ExSampler interfaces to the exsample library.",
     "Events have been sampled using "
     "the ExSample library \\cite{Platzer:2011dr}",
     "%\\cite{Platzer:2011dr}\n"
     "\\bibitem{Platzer:2011dr}\n"
     "S.~Platzer,\n"
     "``ExSample -- A Library for Sampling Sudakov-Type Distributions,''\n"
     "arXiv:1108.6182 [hep-ph].\n"
     "%%CITATION = ARXIV:1108.6182;%%");

  static ThePEG::Parameter<ExSampler,unsigned long> interfacepresampling_points
    ("presampling_points",
     "Set the number of presampling points per cell",
     &ExSampler::presampling_points_, 1000, 0, 0,
     false, false, ThePEG::Interface::lowerlim);

  static ThePEG::Parameter<ExSampler,unsigned long> interfacefreeze_grid
    ("freeze_grid",
     "Set the number of events after which the grid should be frozen",
     &ExSampler::freeze_grid_, 0, 0, 0,
     false, false, ThePEG::Interface::lowerlim);

  static ThePEG::Parameter<ExSampler,double> interfaceefficiency_threshold
    ("efficiency_threshold",
     "Set the efficiency threshold",
     &ExSampler::efficiency_threshold_, .95, 0., 1.,
     false, false, ThePEG::Interface::limited);

  static ThePEG::Parameter<ExSampler,double> interfacegain_threshold
    ("gain_threshold",
     "Set the gain threshold",
     &ExSampler::gain_threshold_, .1, 0., 1.,
     false, false, ThePEG::Interface::limited);

}

