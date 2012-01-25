// -*- C++ -*-

// exsample2/exsampler.cc is part of matchbox
// (C) 2008 Simon Platzer -- sp@particle.uni-karlsruhe.de

#include "exsampler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/ParticleData.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <sstream>

using namespace matchbox;
using namespace exsample;

exsampler::exsampler()
  : presampling_points_(2000),
    freeze_grid_(0),
    efficiency_threshold_(.95),
    gain_threshold_(.1),
    samplers_(), eg_wrappers_(),
    bin_selector_(), integrals_(),
    missing_events_(), oversampling_bins_(),
    sum_weights_(0),
    integral_(0.), integral_abs_(0.), variance_(0.),
    max_integral_(0.), last_bin_(-1),
    last_weight_(0.) {}

exsampler::~exsampler() {}

ThePEG::IBPtr exsampler::clone() const {
  return ThePEG::new_ptr(*this);
}

ThePEG::IBPtr exsampler::fullclone() const {
  return ThePEG::new_ptr(*this);
}

std::string exsampler::process(int bin) const {
  std::ostringstream os("");
  const ThePEG::StandardEventHandler& eh = *eventHandler();
  const ThePEG::StandardXComb& xc = *eh.xCombs()[bin];
  os << xc.mePartonData()[0]->PDGName() << " "
     << xc.mePartonData()[1]->PDGName() << " -> ";
  for ( ThePEG::cPDVector::const_iterator pid =
	  xc.mePartonData().begin() + 2;
	pid != xc.mePartonData().end(); ++pid )
    os << (**pid).PDGName() << " ";
  return os.str();
}

void exsampler::dofinish() {
  ThePEG::SamplerBase::dofinish();

  integral_ = 0.;
  integral_abs_ = 0.;
  variance_ = 0.;
  max_integral_ = 0.;

  int b = 0;

  bool compensating = false;

  for (std::vector<sampler_type>::iterator s = samplers_.begin();
       s != samplers_.end(); ++s, ++b) {

    std::cout << "integrated cross section for "
	      << process(b) << "is ("
	      << samplers_[b].integral() << " +/- "
	      << samplers_[b].integral_uncertainty()
	      << ") nb"
	      << (s->compensating() ? " (*)" : "") << "\n";

    if ( s->compensating() )
      compensating = true;

    s->finalize();
    integral_ += s->integral();
    integral_abs_ += std::abs(s->integral());
    variance_ += s->integral_variance();
    max_integral_ = std::max(max_integral_,s->integral());
  }  

  std::cout << "integrated cross section is ("
	    << integral_ << " +/- "
	    << std::sqrt(variance_) << ") nb\n" << std::flush;

  if ( compensating )
    std::cout << "warning -- samplers marked with (*) are still in compensating mode\n" << std::flush;

}

void exsampler::doinitrun() {
  ThePEG::SamplerBase::doinitrun();
  initialize();
}

void exsampler::persistentOutput(ThePEG::PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << presampling_points_
     << freeze_grid_ << efficiency_threshold_
     << gain_threshold_ << samplers_.size();

  for (std::vector<sampler_type>::const_iterator s = samplers_.begin();
       s != samplers_.end(); ++s)
    s->put(os);

}

void exsampler::persistentInput(ThePEG::PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> presampling_points_
     >> freeze_grid_ >> efficiency_threshold_
     >> gain_threshold_;
  std::size_t nsamplers;
  is >> nsamplers;
  samplers_.resize(nsamplers);

  for (std::vector<sampler_type>::iterator s = samplers_.begin();
       s != samplers_.end(); ++s)
    s->get(is);

}

ThePEG::ClassDescription<exsampler> exsampler::initexsampler;
// Definition of the static class description member.

void exsampler::Init() {

  static ThePEG::ClassDocumentation<exsampler> documentation
    ("There is no documentation for the exsampler class");

  static ThePEG::Parameter<exsampler,unsigned long> interfacepresampling_points
    ("presampling_points",
     "Set the number of presampling points per cell",
     &exsampler::presampling_points_, 2000, 0, 0,
     false, false, ThePEG::Interface::lowerlim);

  static ThePEG::Parameter<exsampler,unsigned long> interfacefreeze_grid
    ("freeze_grid",
     "Set the number of events after which the grid should be frozen",
     &exsampler::freeze_grid_, 0, 0, 0,
     false, false, ThePEG::Interface::lowerlim);

  static ThePEG::Parameter<exsampler,double> interfaceefficiency_threshold
    ("efficiency_threshold",
     "Set the efficiency threshold",
     &exsampler::efficiency_threshold_, .95, 0., 1.,
     false, false, ThePEG::Interface::limited);

  static ThePEG::Parameter<exsampler,double> interfacegain_threshold
    ("gain_threshold",
     "Set the gain threshold",
     &exsampler::gain_threshold_, .1, 0., 1.,
     false, false, ThePEG::Interface::limited);

}

void exsampler::initialize() {

  integral_ = 0.;
  integral_abs_ = 0.;
  variance_ = 0.;
  max_integral_ = 0.;

  if (samplers_.empty())
    samplers_.resize(eventHandler()->nBins());

  eg_wrappers_.resize(eventHandler()->nBins());
  integrals_.resize(eventHandler()->nBins());
  missing_events_.resize(eventHandler()->nBins(),0);

  for (int b = 0; b < eventHandler()->nBins(); ++b) {

    last_bin_ = b;

    eg_wrappers_[b] = eg_exsample2_wrapper(b,eventHandler());
    samplers_[b].function(&eg_wrappers_[b]);

    samplers_[b].sampling_parameters().presampling_points = presampling_points_;
    samplers_[b].sampling_parameters().freeze_grid = freeze_grid_;
    samplers_[b].sampling_parameters().maxtry = eventHandler()->maxLoop();
    samplers_[b].sampling_parameters().efficiency_threshold = efficiency_threshold_;
    samplers_[b].sampling_parameters().gain_threshold = gain_threshold_;

    samplers_[b].initialize(dummy_);
    integrals_[b] = samplers_[b].integral();
    integral_ += samplers_[b].integral();
    integral_abs_ += std::abs(samplers_[b].integral());
    variance_ += samplers_[b].integral_variance();
    max_integral_ = std::max(max_integral_,samplers_[b].integral());

    std::cout << "estimated cross section for " << process(b) << "is ("
	      << samplers_[b].integral() << " +/- "
	      << samplers_[b].integral_uncertainty()
	      << ") nb\n";

  }

  std::cout << "estimated cross section is ("
	    << integral_ << " +/- "
	    << std::sqrt(variance_) << ") nb\n" << std::flush;

  double sum = 0.;

  for (int b = 0; b < eventHandler()->nBins(); ++b) {
    sum += std::abs(samplers_[b].integral());
    bin_selector_.insert(std::make_pair(sum/integral_abs_,b));
  }

}

void exsampler::update() {

  double old_integral_abs = integral_abs_;
  std::vector<double> old_integrals = integrals_;

  integral_ = 0.;
  integral_abs_ = 0.;
  variance_ = 0.;
  max_integral_ = 0.;

  for (std::vector<sampler_type>::iterator s = samplers_.begin();
       s != samplers_.end(); ++s) {

    integral_ += s->integral();
    integral_abs_ += std::abs(s->integral());
    variance_ += s->integral_variance();
    max_integral_ = std::max(max_integral_,s->integral());

  }

  bin_selector_.clear();

  double sum = 0.;

  for (int b = 0; b < eventHandler()->nBins(); ++b) {

    sum += std::abs(samplers_[b].integral());
    integrals_[b] = samplers_[b].integral();
    bin_selector_.insert(std::make_pair(sum/integral_abs_,b));

    long missing = 
      static_cast<long>(round(samplers_[b].stats().accepted()*
			      ((old_integral_abs * std::abs(integrals_[b])) /
			       (integral_abs_ * std::abs(old_integrals[b]))-1.)));

    missing_events_[b] += missing;
    if (missing_events_[b] > 0)
      oversampling_bins_.insert(b);

  }

}

double exsampler::generate() {

  while (true) {

    try {

      if (oversampling_bins_.empty()) {
	last_bin_ = bin_selector_.upper_bound(ThePEG::UseRandom::rnd())->second;
	if (missing_events_[last_bin_] < 0) {
	  ++missing_events_[last_bin_];
	  continue;
	}
      } else {
	last_bin_ = *oversampling_bins_.begin();
	if (--missing_events_[last_bin_] == 0)
	  oversampling_bins_.erase(last_bin_);
      }

      last_weight_ = samplers_[last_bin_].generate(dummy_);
      lastPoint() = samplers_[last_bin_].last_point();
      break;

    } catch (selection_maxtry&) {

      throw maxtry_exception()
	<< "The maximum number of attempts to select a cell was exceeded\n"
	<< "for process " << process(last_bin_)
	<< ThePEG::Exception::eventerror;

    } catch (hit_and_miss_maxtry&) {

      throw maxtry_exception()
	<< "The maximum number of attempts to select an event was exceeded\n"
	<< "process " << process(last_bin_)
	<< ThePEG::Exception::eventerror;

    } catch (generator_update&) {

      update();
      continue;

    } catch (...) {
      throw;
    }

  }

  sum_weights_ += last_weight_;
  return last_weight_;

}

