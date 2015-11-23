// -*- C++ -*-
//
// DipoleSplittingInfo.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleSplittingInfo_H
#define HERWIG_DipoleSplittingInfo_H
//
// This is the declaration of the DipoleIndex and DipoleSplittingInfo classes.
//

#include "ThePEG/PDF/PDF.h"
#include "ThePEG/PDT/ParticleData.h"

#include "Herwig/DipoleShower/Kinematics/DipoleSplittingKinematics.h"

namespace Herwig {

using namespace ThePEG;

class DipoleSplittingKinematics;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 *
 * \brief DipoleIndex is used to index splitting generators
 * for a particular dipole.
 *
 */
class DipoleIndex {

public:

  /**
   * The default constructor.
   */
  DipoleIndex();

  /**
   * The standard constructor
   */
  DipoleIndex(tcPDPtr newEmitter, tcPDPtr newSpectator,
	      const PDF& newEmitterPDF = PDF(), const PDF& newSpectatorPDF = PDF());

public:

  /**
   * Compare for equality.
   */
  bool operator ==(const DipoleIndex& x) const;

  /**
   * Compare for ordering.
   */
  bool operator <(const DipoleIndex& x) const;

  /**
   * Swap emitter and spectator.
   */
  void swap();

  /**
   * Produce a pair of dipole indices given
   * a particle data object for the emission.
   * The ME correction is ignored in the children.
   * The emission is inserted between the emitter
   * and spectator, being a spectator in the first
   * dipole index containing the original emitter,
   * and an emitter in the second dipole, containing
   * the original spectator.
   */
  pair<DipoleIndex,DipoleIndex> split(tcPDPtr) const;

public:

  /**
   * Return the emitter particle data object.
   */
  tcPDPtr emitterData() const { return theEmitterData; }

  /**
   * Return true, if the emitter is an incoming parton
   */
  bool initialStateEmitter() const { return theInitialStateEmitter; }

  /**
   * Return the PDF object associated with the emitter
   */
  const PDF& emitterPDF() const { return theEmitterPDF; }

  /**
   * Return the spectator particle data object.
   */
  tcPDPtr spectatorData() const { return theSpectatorData; }

  /**
   * Return true, if the spectator is an incoming parton
   */
  bool initialStateSpectator() const { return theInitialStateSpectator; }

  /**
   * Return the PDF object associated with the spectator
   */
  const PDF& spectatorPDF() const { return theSpectatorPDF; }

public:

  /**
   * Put information to ostream
   */
  void print(ostream&) const;

private:

  /**
   * The particle data object of the emitter.
   */
  tcPDPtr theEmitterData;

  /**
   * Wether or not the emitter is an incoming parton.
   */
  bool theInitialStateEmitter;

  /**
   * The PDF object for the emitter.
   */
  PDF theEmitterPDF;

  /**
   * The particle data object of the spectator.
   */
  tcPDPtr theSpectatorData;

  /**
   * Wether or not the spectator is an incoming parton.
   */
  bool theInitialStateSpectator;

  /**
   * The PDF object for the spectator.
   */
  PDF theSpectatorPDF;

};

inline ostream& operator << (ostream& os, const DipoleIndex& di) {
  di.print(os);
  return os;
}

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 *
 * \brief DipoleSplittingInfo contains all parameters to generate a full
 * dipole splitting.
 *
 */
class DipoleSplittingInfo {

public:

  /**
   * The default constructor.
   */
  DipoleSplittingInfo();

  /**
   * Destructor
   */
  virtual ~DipoleSplittingInfo() {}

public:

  /**
   * Assign data from another splitting info
   */
  void fill(const DipoleSplittingInfo&);

public:

  /**
   * Return the dipole index
   */
  const DipoleIndex& index() const { return theIndex; }

  /**
   * Return which of the particles
   * in the dipole should be considered emitter (true)
   * and spectator (false)
   */
  const pair<bool,bool>& configuration() const { return theConfiguration; }

  /**
   * Get the configuration marking the spectator
   */
  const pair<bool,bool>& spectatorConfiguration() const { return theSpectatorConfiguration; }

  /**
   * Return the particle data object of the emitter
   * after the splitting.
   */
  tcPDPtr emitterData() const { return theEmitterData; }

  /**
   * Return the particle data object of the emission
   * after the splitting.
   */
  tcPDPtr emissionData() const { return theEmissionData; }

  /**
   * Return the particle data object of the spectator
   * after the splitting.
   */
  tcPDPtr spectatorData() const { return theSpectatorData; }

  /**
   * Return the momentum fraction of the emitter.
   */
  double emitterX() const { return theEmitterX; }

  /**
   * Return the momentum fraction of the spectator.
   */
  double spectatorX() const { return theSpectatorX; }

public:

  /**
   * Return a pointer to the DipoleSplittingKinematics object
   * which is to be used to perform the splitting.
   */
  Ptr<DipoleSplittingKinematics>::tptr splittingKinematics() const { return theSplittingKinematics; }

  /**
   * Return the dipole scale
   */
  Energy scale() const { return theScale; }

  /**
   * Return the pt below which this
   * splitting has been generated.
   */
  Energy hardPt() const { return theHardPt; }

  /**
   * Return the last generated pt
   */
  Energy lastPt() const { return theLastPt; }

  /**
   * Return the last generated momentum fraction.
   */
  double lastZ() const { return theLastZ; }

  /**
   * Return the last generated azimuthal angle.
   */
  double lastPhi() const { return theLastPhi; }

  /**
   * Return the momentum fraction, by which the emitter's
   * momentum fraction should be divided after the splitting.
   */
  double lastEmitterZ() const { return theLastEmitterZ; }

  /**
   * Return the momentum fraction, by which the spectator's
   * momentum fraction should be divided after the splitting.
   */
  double lastSpectatorZ() const { return theLastSpectatorZ; }

  /**
   * Return any additional parameters needed to
   * evaluate the splitting kernel or to generate the 
   * full splitting.
   */
  const vector<double>& lastSplittingParameters() const { return theLastSplittingParameters; }

  /**
   * Return true, if this splitting will terminate
   * the evolution of the dipole considered.
   */
  bool stoppedEvolving() const { return theStoppedEvolving; }

public:

  /**
   * Set the index.
   */
  void index(const DipoleIndex& ind) { theIndex = ind; }

  /**
   * Set the DipoleSplittingKinematics object
   */
  void splittingKinematics(Ptr<DipoleSplittingKinematics>::tptr newSplittingKinematics) {
    theSplittingKinematics = newSplittingKinematics;
  }

  /**
   * Set the particle data object of the emitter
   * after the splitting.
   */
  void emitterData(tcPDPtr p) { theEmitterData = p; }

  /**
   * Set the particle data object of the emission
   * after the splitting.
   */
  void emissionData(tcPDPtr p) { theEmissionData = p; }

  /**
   * Set the particle data object of the spectator
   * after the splitting.
   */
  void spectatorData(tcPDPtr p) { theSpectatorData = p; }

  /**
   * Set the dipole scale
   */
  void scale(Energy s) { theScale = s; }

  /**
   * Set the emitter's momentum fraction
   */
  void emitterX(double x) { theEmitterX = x; }

  /**
   * Set the spectator's momentum fraction
   */
  void spectatorX(double x) { theSpectatorX = x; }

  /**
   * Set the pt below which this
   * splitting has been generated.
   */
  void hardPt(Energy p) { theHardPt = p; }

  /**
   * Set the last generated pt
   */
  void lastPt(Energy p) { theLastPt = p; }

  /**
   * Set the last generated momentum fraction.
   */
  void lastZ(double z) { theLastZ = z; }

  /**
   * Set the last generated azimuthal angle.
   */
  void lastPhi(double p) { theLastPhi = p; }

  /**
   * Set the momentum fraction, by which the emitter's
   * momentum fraction should be divided after the splitting.
   */
  void lastEmitterZ(double z) { theLastEmitterZ = z; }

  /**
   * Set the momentum fraction, by which the spectator's
   * momentum fraction should be divided after the splitting.
   */
  void lastSpectatorZ(double z) { theLastSpectatorZ = z; }

  /**
   * Return the last splitting kernel value encountered.
   */
  double lastValue() const { return theLastValue; }

  /**
   * Set the last splitting kernel value encountered.
   */
  void lastValue(double v) { theLastValue = v; }

  /**
   * Set the last splitting parameters.
   */
  void lastSplittingParameters(const vector<double>& p) { theLastSplittingParameters = p; }

  /**
   * Access the splitting parameters
   */
  vector<double>& splittingParameters() { return theLastSplittingParameters; }

  /**
   * Indicate that this splitting will terminate
   * the evolution of the dipole considered.
   */
  void didStopEvolving() { theStoppedEvolving = true; }

  /**
   * Indicate that this splitting will not terminate
   * the evolution of the dipole considered.
   */
  void continuesEvolving() { theStoppedEvolving = false; }

  /**
   * Reset the configuration.
   */
  void configuration(const pair<bool,bool>& newConfig) { theConfiguration = newConfig; }

  /**
   * Set the configuration marking the spectator
   */
  void spectatorConfiguration(const pair<bool,bool>& conf) { theSpectatorConfiguration = conf; }

public:

  /**
   * Set a pointer to the emitter parton before emission.
   */
  void emitter(tPPtr newEmitter) { theEmitter = newEmitter; }

  /**
   * Set a pointer to the spectator parton before emission.
   */
  void spectator(tPPtr newSpectator) { theSpectator = newSpectator; }

  /**
   * Set a pointer to the emitter parton after emission.
   */
  void splitEmitter(tPPtr newEmitter) { theSplitEmitter = newEmitter; }

  /**
   * Set a pointer to the spectator parton after emission.
   */
  void splitSpectator(tPPtr newSpectator) { theSplitSpectator = newSpectator; }

  /**
   * Set a pointer to the emitted parton.
   */
  void emission(tPPtr newEmission) { theEmission = newEmission; }

  /**
   * Return a pointer to the emitter parton before emission.
   */
  tPPtr emitter() const { return theEmitter; }

  /**
   * Return a pointer to the spectator parton before emission.
   */
  tPPtr spectator() const { return theSpectator; }

  /**
   * Return a pointer to the emitter parton after emission.
   */
  tPPtr splitEmitter() const { return theSplitEmitter; }

  /**
   * Return a pointer to the spectator parton after emission.
   */
  tPPtr splitSpectator() const { return theSplitSpectator; }

  /**
   * Return a pointer to the emitted parton.
   */
  tPPtr emission() const { return theEmission; }

public:

  /**
   * Put information to ostream
   */
  void print(ostream&) const;

private:

  /**
   * The DipoleIndex associated 
   * with this splitting.
   */
  DipoleIndex theIndex;

  /**
   * Flags indicateing which of the particles
   * in the dipole should be considered emitter (true)
   * and spectator (false)
   */
  pair<bool,bool> theConfiguration;

  /**
   * The configuration marking the spectator
   */
  pair<bool,bool> theSpectatorConfiguration;

  /**
   * The particle data object of the emitter
   * after the splitting.
   */
  tcPDPtr theEmitterData;

  /**
   * The particle data object of the emission
   * after the splitting.
   */
  tcPDPtr theEmissionData;

  /**
   * The particle data object of the spectator
   * after the splitting.
   */
  tcPDPtr theSpectatorData;

  /**
   * A pointer to the DipoleSplittingKinematics object
   * which is to be used to perform the splitting.
   */
  Ptr<DipoleSplittingKinematics>::tptr theSplittingKinematics;

  /**
   * The scale for this dipole.
   */
  Energy theScale;

  /**
   * The momentum fraction of the emitter.
   */
  double theEmitterX;

  /**
   * The momentum fraction of the spectator.
   */
  double theSpectatorX;

  /**
   * The pt below which this splitting has 
   * been generated.
   */
  Energy theHardPt;

  /**
   * The last generated pt
   */
  Energy theLastPt;

  /**
   * The last generated momentum fraction.
   */
  double theLastZ;

  /**
   * The last generated azimuthal angle.
   */
  double theLastPhi;

  /**
   * The momentum fraction, by which the emitter's
   * momentum fraction should be divided after the splitting.
   */
  double theLastEmitterZ;

  /**
   * The momentum fraction, by which the spectator's
   * momentum fraction should be divided after the splitting.
   */
  double theLastSpectatorZ;

  /**
   * The last splitting kernel value encountered.
   */
  double theLastValue;

  /**
   * Any additional parameters needed to
   * evaluate the splitting kernel or to generate the 
   * full splitting.
   */
  vector<double> theLastSplittingParameters;

  /**
   * True, if this splitting will terminate
   * the evolution of the dipole considered.
   */
  bool theStoppedEvolving;

  /**
   * A pointer to the emitter parton before emission.
   */
  PPtr theEmitter;

  /**
   * A pointer to the spectator parton before emission.
   */
  PPtr theSpectator;

  /**
   * A pointer to the emitter parton after emission.
   */
  PPtr theSplitEmitter;

  /**
   * A pointer to the spectator parton after emission.
   */
  PPtr theSplitSpectator;

  /**
   * A pointer to the emitted parton.
   */
  PPtr theEmission;

};

inline ostream& operator << (ostream& os, const DipoleSplittingInfo& di) {
  di.print(os);
  return os;
}

}

#endif /* HERWIG_DipoleSplittingInfo_H */
