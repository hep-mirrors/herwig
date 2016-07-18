// -*- C++ -*-
//
// ElectroWeakMatching.h is a part of Herwig - A multi-purpose Monte Carlo event generator
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
//
#ifndef HERWIG_ElectroWeakMatching_H
#define HERWIG_ElectroWeakMatching_H
#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Config/Unitsystem.h"
#include "EWProcess.h"
#include <boost/numeric/ublas/matrix.hpp>

namespace Herwig {
using namespace ThePEG;

namespace ElectroWeakMatching {

  /**
   *  The high energy matching
   */
  boost::numeric::ublas::matrix<Complex>
  electroWeakMatching(Energy mu, Energy2 s, Energy2 t, Energy2 u,
		      Herwig::EWProcess::Process process,
		      bool oneLoop,unsigned int iswap);
}
}

#endif // HERWIG_ElectroWeakMatching_H 
