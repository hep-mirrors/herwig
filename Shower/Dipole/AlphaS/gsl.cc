// -*- C++ -*-

// math/gsl.cc is part of matchbox
// (C) 2008 Simon Platzer -- sp@particle.uni-karlsruhe.de

#include "ThePEG/Utilities/Throw.h"

#include "gsl.h"

using namespace matchbox::gsl;

/// gsl error handler throwing gsl_exception
void error_handler_wrapper (const char * msg,
			    const char *, int, int) {
  ThePEG::Throw<gsl_exception> () << "Matchbox GSL interface : GSL exception : "
				  << msg << ThePEG::Exception::runerror;
}

struct error_handler_resetter_ {
  inline error_handler_resetter_ () {
    gsl_set_error_handler(&error_handler_wrapper);
  }
};

static error_handler_resetter_ error_handler_resetter;

