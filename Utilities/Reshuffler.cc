// -*- C++ -*-
//
// Reshuffler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Reshuffler class.
//

#include <config.h>
#include "Reshuffler.h"
#include "Herwig/Utilities/GSLBisection.h"
#include "Herwig/Shower/ShowerHandler.h"

using namespace Herwig;

void Reshuffler::reshuffle(const PVector& particles,
			   const vector<Energy>& masses) const {

  Energy zero (0.*GeV);
  Lorentz5Momentum Q (zero,zero,zero,zero);
    
  for (PVector::const_iterator p = particles.begin();
       p != particles.end(); ++p) {
    Q += (**p).momentum();
  }

  Boost beta = Q.findBoostToCM();

  vector<Lorentz5Momentum> mbackup;

  bool need_boost = (beta.mag2() > Constants::epsilon);

  if (need_boost) {

    for (PVector::const_iterator p = particles.begin();
	 p != particles.end(); ++p) {
      Lorentz5Momentum mom = (**p).momentum();
      mbackup.push_back(mom);
      (**p).set5Momentum(mom.boost(beta));
    }

  }

  double xi;

  ReshuffleEquation<PVector::const_iterator,vector<Energy>::const_iterator>
    solve (Q.m(),particles.begin(),particles.end(),
	   masses.begin(),masses.end());

  GSLBisection solver(1e-10,1e-8,10000);

  try {
    xi = solver.value(solve,0.0,1.1);
  } catch (GSLBisection::GSLerror) {
    throw Exception("Failed to reshuffle.",Exception::eventerror);
  } catch (GSLBisection::IntervalError) {
    throw Exception("Failed to reshuffle.",Exception::eventerror);
  }

  PVector::const_iterator p = particles.begin();
  vector<Energy>::const_iterator m = masses.begin();
  for (; p != particles.end(); ++p, ++m) {

    Lorentz5Momentum rm = Lorentz5Momentum (xi*(**p).momentum().x(),
					    xi*(**p).momentum().y(),
					    xi*(**p).momentum().z(),
					    sqrt(sqr(*m) +
						 xi*xi*(sqr((**p).momentum().t())-sqr((**p).dataPtr()->mass()))));
     
    rm.rescaleMass();

    if (need_boost) {
      rm.boost(-beta);
    }

    (**p).set5Momentum(rm);

  }

}
