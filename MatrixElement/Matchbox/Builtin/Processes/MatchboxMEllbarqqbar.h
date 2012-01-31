// -*- C++ -*-
//
// MatchboxMEllbarqqbar.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxMEllbarqqbar_H
#define HERWIG_MatchboxMEllbarqqbar_H

#include "ThePEG/Config/Complex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/SpinorHelicity.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/AmplitudeCache.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxMEllbarqqbar implements the
 * generic matrix element for processes with two
 * leptons and two light quarks.
 *
 */
class MatchboxMEllbarqqbar: public SpinorHelicity::AmplitudeCache<int> {

public:

  virtual ~MatchboxMEllbarqqbar() {}

protected:

  /**
   * The Boson mass
   */
  Energy BMass;

  /**
   * The Boson width
   */
  Energy BWidth;

  /**
   * The number of colours.
   */
  double Nc;

  /**
   * The alpha s
   */
  double alphaS;

  /**
   * Lepton QED coupling
   */
  double el;

  /**
   * Lepton vector coupling
   */
  double vl;

  /**
   * Lepton axial vector coupling
   */
  double al;

  /**
   * Quark QED couplings
   */
  pair<double,double> quarkQEDCouplings;

  /**
   * Quark vector couplings
   */
  pair<double,double> quarkVectorCouplings;

  /**
   * Quark axial vector couplings
   */
  pair<double,double> quarkAxialCouplings;

  /**
   * The last calculated inverse photon propagator.
   */
  mutable double Q2;

  /**
   * The last calculated inverse Z propagator.
   */
  mutable Complex Q2BW;

  /**
   * The last quark QED coupling
   */
  mutable double eq;

  /**
   * The last quark vector coupling
   */
  mutable double vq;

  /**
   * The last quark axial vector coupling
   */
  mutable double aq;

  /**
   * Initialize couplings and masses. This default version
   * assumes Z/gamma exchange.
   */
  virtual void doinit(const HandlerBase& h);

  /**
   * Write out persistently
   */
  void persistentOutput(PersistentOStream & os) const;
						      
  /**
   * Read in persistently
   */
  void persistentInput(PersistentIStream & is);

  /**
   * Prepare for given momenta.
   */
  void prepare(const Lorentz5Momentum& pl, const Lorentz5Momentum& plbar,
	       const Lorentz5Momentum& pq, const Lorentz5Momentum& pqbar,
	       Energy2 sHat, 
	       cPDPtr lData, cPDPtr lbarData, 
	       cPDPtr qData, cPDPtr qbarData) const;

  /**
   * Prepare for given momenta.
   */
  void prepare(const Lorentz5Momentum& pl, const Lorentz5Momentum& plbar,
	       const Lorentz5Momentum& pq, const Lorentz5Momentum& pqbar,
	       Energy2 sHat, 
	       cPDPtr lData, cPDPtr qData) const {
    prepare(pl,plbar,pq,pqbar,sHat,lData,lData,qData,qData);
  }

  /**
   * Return the amplitude squared summed over helicities
   */
  double evaluateME2(bool photon = true) const;

};

}

#endif // HERWIG_MatchboxMEllbarqqbar_H
