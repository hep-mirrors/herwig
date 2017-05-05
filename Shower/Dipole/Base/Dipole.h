// -*- C++ -*-
//
// Dipole.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Dipole_H
#define HERWIG_Dipole_H
//
// This is the declaration of the Dipole class.
//

#include "Herwig/Shower/Dipole/Kinematics/DipoleSplittingKinematics.h"
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer, Stephen Webster, Johannes Bellm
 *
 * \brief The Dipole class is used by the dipole shower to
 * represent a dipole of two coloured partons.
 *
 */
class Dipole {

public:

  /**
   * The default constructor
   */
  Dipole();
  
  /**
   * The standard constructor
   */
  Dipole(const pair<PPtr,PPtr>& newParticles,
	 const pair<PDF,PDF>& newPDFs,
	 pair<double,double> newFractions,
	 pair<Energy,Energy> newScales);

  /**
   * The standard constructor
   */
  Dipole(const pair<PPtr,PPtr>& newParticles,
	 const pair<PDF,PDF>& newPDFs,
	 pair<double,double> newFractions,
	 const pair<bool,bool> decaying = pair<bool,bool>(false,false),
	 pair<Energy,Energy> newScales = pair<Energy,Energy>(ZERO,ZERO));


public:

  /**
   * Get the left particle.
   */
  tPPtr leftParticle() const { return theParticles.first; }

  /**
   * Get the right particle.
   */
  tPPtr rightParticle() const { return theParticles.second; }

  /**
   * Get the left PDF.
   */
  const PDF& leftPDF() const { return thePDFs.first; }

  /**
   * Get the right PDF.
   */
  const PDF& rightPDF() const { return thePDFs.second; }

  /**
   * Get the left fraction.
   */
  double leftFraction() const { return theFractions.first; }

  /**
   * Get the right fraction.
   */
  double rightFraction() const { return theFractions.second; }

  /**
   * Get the bool indicating
   * incoming decay for the left
   * particle, for debugging only.
   */
  bool leftDecaying() { return theDecaying.first; }

  /**
   * Get the bool indicating
   * incoming decay for the right
   * particle, for debugging only.
   */
  bool rightDecaying() { return theDecaying.second; }


  /**
   * Set the left particle.
   */
  void leftParticle(PPtr p) { theParticles.first = p; }

  /**
   * Set the right particle.
   */
  void rightParticle(PPtr p) { theParticles.second = p; }

  /**
   * Set the left PDF
   */
  void leftPDF(const PDF& p) { thePDFs.first = p; }

  /**
   * Set the right PDF
   */
  void rightPDF(const PDF& p) { thePDFs.second = p; }

  /**
   * Set the momentum fraction for the left particle.
   */
  void leftFraction(double x) { theFractions.first = x; }

  /**
   * Set the momentum fraction for the right particle.
   */
  void rightFraction(double x) { theFractions.second = x; }

  /**
   * Get the scale for the left particle.
   */
  Energy leftScale() const { return theScales.first; }

  /**
   * Set the scale for the left particle.
   */
  void leftScale(Energy s) { theScales.first = s; }

  /**
   * Get the scale for the right particle.
   */
  Energy rightScale() const { return theScales.second; }

  /**
   * Set the scale for the right particle.
   */
  void rightScale(Energy s) { theScales.second = s; }

  /**
   * Set the decayed particle indicator
   * for the left particle
   */
  void leftDecaying(bool decaying) { theDecaying.first = decaying; }

  /**
   * Set the decayed particle indicator
   * for the right particle
   */
  void rightDecaying(bool decaying) { theDecaying.second = decaying; }

  /**
   * Update information, if modified.
   */
  void update();

public:

  /**
   * Return the dipole index for the selected
   * emitter-spectator assignment.
   */
  const DipoleIndex& index(pair<bool,bool> conf) const {
    return conf.first ? theIndices.first : theIndices.second;
  }

  /**
   * Return the emitter particle for the
   * selected configuration.
   */
  tPPtr emitter(pair<bool,bool> conf) const {
    return conf.first ? theParticles.first : theParticles.second;
  }

  /**
   * Return the spectator particle for the
   * selected configuration.
   */
  tPPtr spectator(pair<bool,bool> conf) const {
    return conf.first ? theParticles.second : theParticles.first;
  }

  /**
   * Return the scale associated to the emitter
   * for the selected configuration.
   */
  Energy emitterScale(pair<bool,bool> conf) const {
    return conf.first ? theScales.first : theScales.second;
  }

  /**
   * Set the scale associated to the emitter
   * for the selected configuration.
   */
  void emitterScale(pair<bool,bool> conf, Energy scale) {
    (conf.first ? theScales.first : theScales.second) = scale;
  }

  /**
   * Return the momentum fraction of the emitter
   * for the selected configuration.
   */
  double emitterX(pair<bool,bool> conf) const {
    return conf.first ? theFractions.first : theFractions.second;
  }

  /**
   * Return the PDF of the emitter
   * for the selected configuration.
   */
  const PDF& emitterPDF(pair<bool,bool> conf) const {
    return conf.first ? thePDFs.first : thePDFs.second;
  }

  /**
   * Return the momentum fraction of the spectator
   * for the selected configuration.
   */
  double spectatorX(pair<bool,bool> conf) const {
    return conf.first ? theFractions.second : theFractions.first;
  }

  /**
   * Return the PDF of the spectator
   * for the selected configuration.
   */
  const PDF& spectatorPDF(pair<bool,bool> conf) const {
    return conf.first ? thePDFs.second : thePDFs.first;
  }

public:

  /**
   * Split this dipole according to the given splitting.
   * If colourSpectator is true, do not change the spectator.
   */
  pair<Dipole,Dipole> split (DipoleSplittingInfo& dsplit,
			     bool colourSpectator) const;

  /**
   * As split, but without touching the event record.
   * Needed to produce a phase space point as it would 
   * be after calling split.
   */
  void tmpsplit (DipoleSplittingInfo& dsplit,
                             bool colourSpectator) const;


  /**
   * Produce a new spectator according to the
   * given splitting.
   */
  void recoil (DipoleSplittingInfo& dsplit);

public:

  /**
   * Put information to ostream
   */
  void print(ostream&) const;

private:

  /**
   * The particles forming the dipole
   */
  pair<PPtr,PPtr> theParticles;

  /**
   * The PDF objects.
   */
  pair<PDF,PDF> thePDFs;

  /**
   * The momentum fractions associated
   * to the incoming particles
   */
  pair<double,double> theFractions;

  /**
   * The dipole indices, if the first or second particle
   * is considered as emitter.
   */
  pair<DipoleIndex,DipoleIndex> theIndices;

  /**
   * The scale associated to the first and second
   * particle, respectively.
   */
  pair<Energy,Energy> theScales;

  /**
   * Indicates if either the first or the second parton 
   * is incoming to a decay.
   */
  pair<bool,bool> theDecaying;

};

inline ostream& operator << (ostream& os, const Dipole& di) {
  di.print(os);
  return os;
}

}

#endif /* HERWIG_Dipole_H */
