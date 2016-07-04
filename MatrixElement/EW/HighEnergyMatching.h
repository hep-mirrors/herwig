// -*- C++ -*-
//
// HighEnergyMatching.h is a part of Herwig - A multi-purpose Monte Carlo event generator
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
//
#ifndef HERWIG_HighEnergyMatching_H
#define HERWIG_HighEnergyMatching_H
#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Config/Unitsystem.h"
#include "EWProcess.h"
#include <boost/numeric/ublas/matrix.hpp>

namespace Herwig {
using namespace ThePEG;

namespace HighEnergyMatching {

  /**
   *  The high energy matching
   */
  boost::numeric::ublas::matrix<complex<InvEnergy2> > 
    highEnergyMatching(Energy highScale, 
		       Energy2 s, Energy2 t, Energy2 u,
		       Herwig::EWProcess::Process process,
		       bool oneLoop, bool includeAlphaS2);
  
  /**
   *  Spin\f$\frac12\f$
   */
  boost::numeric::ublas::matrix<complex<InvEnergy2> > 
    SpinHalfMatching(Energy highScale, 
		     Energy2 s, Energy2 t, Energy2 u,
		     EWProcess::Process process,
		     bool oneLoop, bool includeAlphaS2);
  
  /**
   *  Spin\f$1\f$
   */
  boost::numeric::ublas::matrix<complex<InvEnergy2> > 
    Spin1HighMatching(Energy highScale, 
		      Energy2 s, Energy2 t, Energy2 u,
		      EWProcess::Process process,
		      bool oneLoop, bool includeAlphaS2);
  /**
   *  Spin\f$0\f$
   */
  boost::numeric::ublas::matrix<complex<InvEnergy2> > 
    Spin0HighMatching(Energy highScale, 
		      Energy2 s, Energy2 t, Energy2 u,
		      EWProcess::Process process,
		      bool oneLoop, bool includeAlphaS2);
}
}

#endif // HERWIG_HighEnergyMatching_H 
