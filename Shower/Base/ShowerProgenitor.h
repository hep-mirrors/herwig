// -*- C++ -*-
//
// ShowerProgenitor.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerProgenitor_H
#define HERWIG_ShowerProgenitor_H
//
// This is the declaration of the ShowerProgenitor struct.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Shower/ShowerConfig.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ShowerProgenitor.fh"
#include "ThePEG/PDF/BeamParticleData.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  A struct to store information on the perturbative particle which 
 *  initiates a shower
 */
class ShowerProgenitor : public Base {

/**
 *  Typedef for the BeamParticleData objects
 */
typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

public:

  /**
   *  Constructor for the class
   * @param original The original particle
   * @param copy The colour isolated copy
   * @param particle The ShowerPArticle copy
   * @param pT The \f$p_t\f$ of the hardest emission
   * @param emitted Whether or not the particle has radiated
   */
  ShowerProgenitor(PPtr original,PPtr copy, ShowerParticlePtr particle,
		   Energy pT=ZERO,bool emitted=false)
    : _original(original), _copy(copy), _perturbative(true),
      _particle(particle), _highestpT(pT), _maxpT(Constants::MaxEnergy), 
      _maxHardPt(ZERO), _hasEmitted(emitted) {
    // get the BeamParticleData object
    if ( original->parents().empty() ) {
      _beam=dynamic_ptr_cast<tcBeamPtr>(original->dataPtr());
    } 
    else {
      _beam=dynamic_ptr_cast<tcBeamPtr>(original->parents()[0]->dataPtr());
    }
  }
  
  /**
   *  Access to the particle
   */
  ShowerParticlePtr progenitor() const { return _particle; }

  /**
   *  Set the particle
   */
  void progenitor(ShowerParticlePtr in) { _particle=in; }

  /**
   *  Access to the original particle
   */
  PPtr original() const { return _original; }

  /**
   *  Access to the colour isolated copy
   */
  PPtr copy() const { return _copy; }

  /**
   * Set the copy
   */
  void copy(PPtr in) { _copy=in; }

  /**
   *  Whether the particle came from the hard process or was added by
   *  the matrix element correction
   */
  bool perturbative() const { return _perturbative; }

  /**
   *  Whether the particle came from the hard process or was added by
   *  the matrix element correction
   */
  void perturbative(bool in) { _perturbative=in; }

  /**
   *  Set/Get methods for the hardest \f$p_T\f$ so far
   */
  //@{
  /**
   *  Access the \f$p_T\f$ of the hardest emission so far
   */
  Energy highestpT() const { return _highestpT; }

  /**
   *  Set the \f$p_T\f$ of the hardest emission so far
   */
  void highestpT(Energy in) { _highestpT=in; }
  //@}

  /**
   *  Set/Get methods for the maximum \f$p_T\f$ 
   */
  //@{
  /**
   *  Access the maximum \f$p_T\f$ for radiation
   */
  Energy maximumpT() const { return _maxpT; }

  /**
   *  Set the maximum \f$p_T\f$ for radiation
   */
  void maximumpT(Energy in) { _maxpT=in; }
  //@}

  /**
   * Set/Get methods for whether the particle has radiated
   */
  //@{
  /**
   *  Access the maximum hard \f$p_T\f$, given by the hard process
   */
  Energy maxHardPt() const { return _maxHardPt; }

  /**
   *  Set the maximum hard \f$p_T\f$, given by the hard process
   */
  void maxHardPt(Energy in) { _maxHardPt = in; }

  /**
   *  Has this particle radiated
   */
  bool hasEmitted() const { return _hasEmitted; }

  /**
   *  Set whether or not this particle has radiated
   */
  void hasEmitted(bool in) { _hasEmitted=in; }
  //@}

  /**
   *  The id of the particle
   */
  long id() const { return _particle->id(); }

  /**
   *  The BeamParticleData object
   */
  tcBeamPtr beam() { return _beam; }

private:

  /**
   *  Pointer to the original particle
   */
  PPtr _original;

  /**
   *  Pointer to the colour isolated copy of the original particle
   */
  PPtr _copy;

  /**
   *  Whether the particle came from the hard process or was added by
   *  the matrix element correction
   */
  bool _perturbative;

  /**
   *  Pointer to the ShowerParticle
   */
  ShowerParticlePtr _particle;

  /**
   *  Highest \f$p_T\f$ emitted in the shower from this particle
   */
  Energy _highestpT;

  /**
   *  Maximum allowed \f$p_T\f$ for emission from this particle
   */
  Energy _maxpT;

  /**
   *  maximum hard \f$p_T\f$ from the hard process
   */
  Energy _maxHardPt;

  /**
   *  Has there been radiation
   */
  bool _hasEmitted;

  /**
   *  The BeamParticleData object
   */
  tcBeamPtr _beam;

};
}

#endif /* HERWIG_ShowerProgenitor_H */
