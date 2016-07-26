// -*- C++ -*-

// math/gsl.h is part of matchbox
// (C) 2008 Simon Platzer -- sp@particle.uni-karlsruhe.de

#ifndef matchbox_math_gsl_h
#define matchbox_math_gsl_h

#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"

#include "ThePEG/Utilities/Exception.h"

namespace matchbox {

  namespace gsl {

    /// exception class for GSL problems      
    struct gsl_exception : public ThePEG::Exception { };

    /// wrapper araound the bisection root solver
    template<class Function, unsigned long MaxIterations>
    struct bisection_root_solver {

      /// constructor -- allocates the solver
      inline bisection_root_solver (const Function& thef) : f(thef) {
	s = gsl_root_fsolver_alloc (gsl_root_fsolver_bisection);
      }

      /// destructor -- frees the solver
      ~bisection_root_solver () {
	gsl_root_fsolver_free (s);
      }

      /// solve for root given initial interval
      double solve (std::pair<double,double> interval, double precision = .000001);

      /// function object representing the equation to be solved
      Function f;

      /// the gsl solver used
      gsl_root_fsolver * s;

    };

  }

}

#include "gsl.tcc"

#endif // matchbox_math_gsl_h
