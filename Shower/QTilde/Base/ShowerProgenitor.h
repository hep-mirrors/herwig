// -*- C++ -*-
//
// ShowerProgenitor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerProgenitor_H
#define HERWIG_ShowerProgenitor_H
//
// This is the declaration of the ShowerProgenitor struct.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig/Shower/QTilde/ShowerConfig.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"
#include "ShowerProgenitor.fh"
#include "ThePEG/PDF/BeamParticleData.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  A struct to store information on the perturbative particle which 
 *  initiates a shower
 */
class ShowerProgenitor : public Base {

public:

  /**
   *  Enum for the reconstruction state of this progentitor
   */
  enum Reconstructed { notReconstructed=0, done, dontReconstruct};


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
      _particle(particle), _highestpT(pT), 
      _maxHardPt(ZERO), _hardScale(ZERO), _hasEmitted(emitted),
      _reconstructed(notReconstructed) {
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
  Energy maximumpT(ShowerInteraction type) const {
    assert(type!=ShowerInteraction::Both);
    map<ShowerInteraction,Energy>::const_iterator it = _maxpT.find(type);
    return it !=_maxpT.end() ? it->second : Constants::MaxEnergy; 
  }

  /**
   *  Set the maximum \f$p_T\f$ for radiation
   */
  void maximumpT(Energy in,ShowerInteraction type) {
    _maxpT[type]=in; }
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
   *  Access the relevant hard scale to be used in profile scales; usually
   *  this is taken to be the maximum pt, except for other choices such as
   *  hfact.
   */
  Energy hardScale() const { return _hardScale; }

  /**
   *  Set the relevant hard scale to be used in profile scales; usually
   *  this is taken to be the maximum pt, except for other choices such as
   *  hfact.
   */
  void hardScale(Energy in) { _hardScale = in; }

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

  /**
   *  Whether or not the recon has been performed
   */
  Reconstructed reconstructed() const {return _reconstructed;}

  /**
   *  Whether or not the recon has been performed
   */
  void reconstructed(Reconstructed recon) {_reconstructed = recon;}

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
  map<ShowerInteraction,Energy> _maxpT;

  /**
   *  maximum hard \f$p_T\f$ from the hard process
   */
  Energy _maxHardPt;

  /**
   *  The relevant hard scale to be used in profile scales; usually
   *  this is taken to be the maximum pt, except for other choices such as
   *  hfact.
   */
  Energy _hardScale;

  /**
   *  Has there been radiation
   */
  bool _hasEmitted;

  /**
   *  The BeamParticleData object
   */
  tcBeamPtr _beam;

  /**
   *  Whether or not the reconstruction has been performed
   */
  Reconstructed _reconstructed;

};
}

#endif /* HERWIG_ShowerProgenitor_H */
