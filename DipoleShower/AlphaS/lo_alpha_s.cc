// -*- C++ -*-

// couplings/lo_alpha_s.cc is part of matchbox
// (C) 2008 Simon Platzer -- sp@particle.uni-karlsruhe.de

#include "lo_alpha_s.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace matchbox;

lo_alpha_s::lo_alpha_s()
  : alpha_s(), freezing_scale_(1.*GeV) {}

lo_alpha_s::~lo_alpha_s() {}

IBPtr lo_alpha_s::clone() const {
  return new_ptr(*this);
}

IBPtr lo_alpha_s::fullclone() const {
  return new_ptr(*this);
}

void lo_alpha_s::persistentOutput(PersistentOStream & os) const {
  os << ounit(freezing_scale_,GeV);
}

void lo_alpha_s::persistentInput(PersistentIStream & is, int) {
  is >> iunit(freezing_scale_,GeV);
}

ClassDescription<lo_alpha_s> lo_alpha_s::initlo_alpha_s;
// Definition of the static class description member.

void lo_alpha_s::Init() {

  static ClassDocumentation<lo_alpha_s> documentation
    ("LO running alpha_s");


  static Parameter<lo_alpha_s,Energy> interfacefreezing_scale
    ("freezing_scale",
     "Freeze alpha_s below given scale",
     &lo_alpha_s::freezing_scale_, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

}

double lo_alpha_s::operator () (Energy2 scale,
				Energy2 lambda2,
				unsigned int nf) const {

  if (scale < sqr(freezing_scale_)) {
    scale = sqr(freezing_scale_);
    nf = active_flavours(scale);
    lambda2 = lambda_squared(nf);
  }

  double beta0 = (33.-2.*nf)/(12.*Constants::pi);
	
  return 1./(beta0*log(scale/lambda2));

}
