// -*- C++ -*-

// couplings/nlo_alpha_s.cc is part of matchbox
// (C) 2008 Simon Platzer -- sp@particle.uni-karlsruhe.de

#include "nlo_alpha_s.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace matchbox;

nlo_alpha_s::nlo_alpha_s()
  : alpha_s(), freezing_scale_(1.*GeV), exact_evaluation_(true) {
}

nlo_alpha_s::~nlo_alpha_s() {}

IBPtr nlo_alpha_s::clone() const {
  return new_ptr(*this);
}

IBPtr nlo_alpha_s::fullclone() const {
  return new_ptr(*this);
}

void nlo_alpha_s::persistentOutput(PersistentOStream & os) const {
  os << ounit(freezing_scale_,GeV) << exact_evaluation_;
}

void nlo_alpha_s::persistentInput(PersistentIStream & is, int) {
  is >> iunit(freezing_scale_,GeV) >> exact_evaluation_;
}

ClassDescription<nlo_alpha_s> nlo_alpha_s::initnlo_alpha_s;
// Definition of the static class description member.

void nlo_alpha_s::Init() {

  static ClassDocumentation<nlo_alpha_s> documentation
    ("NLO running alpha_s");


  static Parameter<nlo_alpha_s,Energy> interfacefreezing_scale
    ("freezing_scale",
     "Freeze alpha_s below given scale",
     &nlo_alpha_s::freezing_scale_, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);


  static Switch<nlo_alpha_s,bool> interfaceexact_evaluation
    ("exact_evaluation",
     "Wether to exactly evaluate the running or use running for large scales",
     &nlo_alpha_s::exact_evaluation_, true, true, false);
  static SwitchOption interfaceexact_evaluationexact
    (interfaceexact_evaluation,
     "exact",
     "Perform exact evaluation",
     true);
  static SwitchOption interfaceexact_evaluationlarge_scale
    (interfaceexact_evaluation,
     "large_scale",
     "Perform approximate evaluation for large scales",
     false);

}

double nlo_alpha_s::operator () (Energy2 scale,
				 Energy2 lambda2,
				 unsigned int nf) const {

  if (scale < sqr(freezing_scale_)) {
    scale = sqr(freezing_scale_);
    nf = active_flavours(scale);
    lambda2 = lambda_squared(nf);
  }

  double beta0 = (33.-2.*nf)/(12.*Constants::pi);
  double beta1 = (153.-19.*nf)/(24.*sqr(Constants::pi));

  if (exact_evaluation_) {

    rg_solver().f.slog = log(scale/lambda2);
    rg_solver().f.nf = nf;

    double slog = rg_solver().f.slog;

    double center = 
      (1./(beta0*slog))*
      (1. - (beta1/sqr(beta0)) * log(slog)/slog +
       sqr(beta1/(sqr(beta0)*slog)) * (sqr(log(slog)-.5) - 5./4.));

    return rg_solver().solve(make_pair(.5*center,1.5*center));

  } else {

    double slog = log(scale/lambda2);
	
    return 
      (1./(beta0*slog))*
      (1. - (beta1/sqr(beta0)) * log(slog)/slog +
       sqr(beta1/(sqr(beta0)*slog)) * (sqr(log(slog)-.5) - 5./4.));

  }

  return 0.;

}
