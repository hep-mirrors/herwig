// -*- C++ -*-

// couplings/alpha_s.cc is part of matchbox
// (C) 2008 Simon Platzer -- sp@particle.uni-karlsruhe.de

#include "alpha_s.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/Repository.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace matchbox;

alpha_s::alpha_s()
  : AlphaSBase(), min_active_flavours_(3), max_active_flavours_(6),
    matched_(false), scale_factor_(1.), quark_masses_squared_(),
    lambda_squared_(), alpha_s_in_(.1176), scale_in_(91.1876*GeV),
    lambda_range_(1.*MeV2,1.e6*MeV2), fixed_(false) {
}

alpha_s::~alpha_s() {}

void alpha_s::persistentOutput(PersistentOStream & os) const {
  os << min_active_flavours_ << max_active_flavours_ << matched_ << scale_factor_;
  for (size_t f = 0; f < 7; ++f)
    os << ounit(quark_masses_squared_[f],MeV2)
       << ounit(lambda_squared_[f],MeV2);
  for (size_t f = 0; f < 6; ++f)
    os << ounit(nfvector[f],MeV2);
  os << alpha_s_in_ << ounit(scale_in_,GeV)
     << ounit(lambda_range_.first,MeV2) << ounit(lambda_range_.second,MeV2)
     << fixed_;
}

void alpha_s::persistentInput(PersistentIStream & is, int) {
  is >> min_active_flavours_ >> max_active_flavours_ >> matched_ >> scale_factor_;
  for (size_t f = 0; f < 7; ++f)
    is >> iunit(quark_masses_squared_[f],MeV2)
       >> iunit(lambda_squared_[f],MeV2);
  for (size_t f = 0; f < 6; ++f)
     is >> iunit(nfvector[f],MeV2);
  is >> alpha_s_in_ >> iunit(scale_in_,GeV)
     >> iunit(lambda_range_.first,MeV2) >> iunit(lambda_range_.second,MeV2)
     >> fixed_;
}

AbstractClassDescription<alpha_s> alpha_s::initalpha_s;
// Definition of the static class description member.

void alpha_s::Init() {

  static ClassDocumentation<alpha_s> documentation
    ("Base class for strong coupoling as used in matchbox");


  static Parameter<alpha_s,unsigned int> interfacemin_active_flavours
    ("min_active_flavours",
     "Minimum number of active flavours",
     &alpha_s::min_active_flavours_, 3, 0, 6,
     true, false, Interface::limited);

  static Parameter<alpha_s,unsigned int> interfacemax_active_flavours
    ("max_active_flavours",
     "Maximum number of active flavours",
     &alpha_s::max_active_flavours_, 6, 0, 6,
     true, false, Interface::limited);


  static Parameter<alpha_s,double> interfaceinput_alpha_s
    ("input_alpha_s",
     "alpha_s value at input scale",
     &alpha_s::alpha_s_in_, .1176, 0.0, 1.0,
     true, false, Interface::limited);


  static Parameter<alpha_s,Energy> interfaceinput_scale
    ("input_scale",
     "Input scale for alpha_s value",
     &alpha_s::scale_in_, GeV, 91.1876*GeV, 0.*GeV, 0.*GeV,
     true, false, Interface::lowerlim);


  static Command<alpha_s> interfacecheck
    ("check",
     "check",
     &alpha_s::check, false);

  static Parameter<alpha_s,double> interfacescale_factor
    ("scale_factor",
     "scale factor for argument",
     &alpha_s::scale_factor_, 1., 0.0, 100.0,
     true, false, Interface::limited);


  static Switch<alpha_s,bool> interfacefixed
    ("fixed",
     "",
     &alpha_s::fixed_, false, false, false);
  static SwitchOption interfacefixedYes
    (interfacefixed,
     "Yes",
     "",
     true);
  static SwitchOption interfacefixedNo
    (interfacefixed,
     "No",
     "",
     false);
  
}

string alpha_s::check (string args) {

  istringstream argin(args);

  double Q_low, Q_high;
  long n_steps;

  argin >> Q_low >> Q_high >> n_steps;

  string fname;
  argin >> fname;

  generator()->log() <<  "checking alpha_s in range [" << Q_low << "," << Q_high << "] GeV in "
	    << n_steps << " steps.\nResults are written to " << fname << "\n";

  double step_width = (Q_high-Q_low)/n_steps;

  match_thresholds();

  generator()->log() <<  "threshold matching results:\n"
	    << "(m_Q^2 -> Lambda^2) / GeV^2 for dynamic flavours in range ["
	    << min_active_flavours_ << "," << max_active_flavours_ << "]\n";

  for (size_t f = 0; f < 7; ++f) {
    generator()->log() <<  (quark_masses_squared_[f]/GeV2) << " " 
	      << (lambda_squared_[f]/GeV2) << "\n";
  }

  ofstream out (fname.c_str());

  for (long k = 0; k <= n_steps; ++k) {

    Energy Q = Q_low*GeV + k*step_width*GeV;

    out << (Q/GeV) << " " << (operator () (Q*Q)) << "\n";

  }

  return "alpha_s check finished";

}

void alpha_s::match_thresholds () {

  if (matched_)
    return;

  // get the quark masses
  quark_masses_squared_[0] = 0.*MeV2;

  for (long f = 1; f < 7; ++f) {
    if ( quarkMasses().empty() )
      quark_masses_squared_[static_cast<size_t>(f)]
	= sqr(getParticleData(f)->mass());
    else
      quark_masses_squared_[static_cast<size_t>(f)]
	= sqr(quarkMasses()[static_cast<size_t>(f-1)]);
  }

  if ( quark_masses_squared_[1] > quark_masses_squared_[2] )
    swap(quark_masses_squared_[1],quark_masses_squared_[2]);

  unsigned int active_at_input = active_flavours(sqr(scale_in_));

  // solve for input lambda
  solve_input_lambda<alpha_s> input_equation (this,active_at_input,alpha_s_in_,sqr(scale_in_));

  gsl::bisection_root_solver<solve_input_lambda<alpha_s>,100> input_solver(input_equation);

  lambda_squared_[active_at_input] =
    MeV2 *
  input_solver.solve({lambda_range_.first/MeV2,lambda_range_.second/MeV2});

  // get lambdas down to min active flavours
  unsigned int below = active_at_input;

  while (below > min_active_flavours_) {

    solve_lambda_below<alpha_s> match_equation (this,below,
						lambda_squared_[below],
						quark_masses_squared_[below]);
    gsl::bisection_root_solver<solve_lambda_below<alpha_s>,100> match_solver(match_equation);
    lambda_squared_[below-1] =
      MeV2 *
    match_solver.solve({lambda_range_.first/MeV2,lambda_range_.second/MeV2});

    --below;
  }

  // get lambdas up to max active flavours
  unsigned int above = active_at_input;

  while (above < max_active_flavours_) {
    solve_lambda_above<alpha_s> match_equation (this,above,
						lambda_squared_[above],
						quark_masses_squared_[above+1]);
    gsl::bisection_root_solver<solve_lambda_above<alpha_s>,100> match_solver(match_equation);
    lambda_squared_[above+1] =
      MeV2 *match_solver.solve({lambda_range_.first/MeV2,lambda_range_.second/MeV2});
    ++above;
  }

  if (min_active_flavours_ > 0) {
    for (size_t f = 0; f < min_active_flavours_; ++f) {
      lambda_squared_[f] = lambda_squared_[min_active_flavours_];
    }
  }

  if (max_active_flavours_ < 6) {
    for (size_t f = max_active_flavours_+1; f < 7; ++f) {
      lambda_squared_[f] = lambda_squared_[max_active_flavours_];
    }
  }

  matched_ = true;

  return;

}
