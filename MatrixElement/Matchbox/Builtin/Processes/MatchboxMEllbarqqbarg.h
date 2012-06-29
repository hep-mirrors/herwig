// -*- C++ -*-
//
// MatchboxMEllbarqqbarg.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxMEllbarqqbarg_H
#define HERWIG_MatchboxMEllbarqqbarg_H

#include "Herwig++/MatrixElement/Matchbox/Builtin/Processes/MatchboxMEllbarqqbar.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxMEllbarqqbarg implements the
 * generic matrix element for processes with two
 * leptons, two light quarks and a gluon.
 *
 */
class MatchboxMEllbarqqbarg: public MatchboxMEllbarqqbar {

protected:

  /**
   * Initialize couplings and masses. This default version
   * assumes Z/gamma exchange.
   */
  virtual void doinit(const HandlerBase& h) {
    MatchboxMEllbarqqbar::doinit(h);
    nPoints(5);
  }

  /**
   * Write out persistently
   */
  void persistentOutput(PersistentOStream & os) const {
    MatchboxMEllbarqqbar::persistentOutput(os);
  }
						      
  /**
   * Read in persistently
   */
  void persistentInput(PersistentIStream & is) {
    MatchboxMEllbarqqbar::persistentInput(is);
    nPoints(5);
  }

  /**
   * Prepare for given momenta.
   */
  void prepare(const Lorentz5Momentum& pl, const Lorentz5Momentum& plbar,
	       const Lorentz5Momentum& pq, const Lorentz5Momentum& pqbar,
	       const Lorentz5Momentum& pg,
	       Energy2 sHat, 
	       cPDPtr lData, cPDPtr lbarData,
	       cPDPtr qData, cPDPtr qbarData) const {
    MatchboxMEllbarqqbar::prepare(pl,plbar,pq,pqbar,sHat,lData,lbarData,qData,qbarData);
    momentum(4,pg,false);
  }

  /**
   * Prepare for given momenta.
   */
  void prepare(const Lorentz5Momentum& pl, const Lorentz5Momentum& plbar,
	       const Lorentz5Momentum& pq, const Lorentz5Momentum& pqbar,
	       const Lorentz5Momentum& pg,
	       Energy2 sHat, cPDPtr lData, cPDPtr qData) const {
    prepare(pl,plbar,pq,pqbar,pg,sHat,lData,lData,qData,qData);
  }

  /**
   * Return the amplitude squared summed over helicities
   */
  double evaluateME2(bool photon = true) const;

};

}

#endif // HERWIG_MatchboxMEllbarqqbarg_H
