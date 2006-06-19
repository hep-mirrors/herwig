// -*- C++ -*-
#ifndef HERWIG_ShowerProgenitor_H
#define HERWIG_ShowerProgenitor_H
//
// This is the declaration of the ShowerProgenitor struct.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Shower/Kinematics/ShowerParticle.h"
#include "ShowerConfig.h"

#include "ShowerProgenitor.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  A struct to store information on the perturbative particle which 
 *  initiates a shower
 */
class ShowerProgenitor : public Base {

public:

  /**
   *  Constructor for the class
   * @param original The original particle
   * @param copy The colour isolated copy
   * @param particle The ShowerPArticle copy
   * @param pT The \f$p_t\f$ of the hardest emission
   * @param emitted Whether or not the particle has radiated
   */
  inline ShowerProgenitor(PPtr original,PPtr copy, ShowerParticlePtr particle,
			  Energy pT=0.,bool emitted=false);

  /**
   *  Access to the particle
   */
  inline ShowerParticlePtr progenitor() const;

  /**
   *  Set the particle
   */
  inline void progenitor(ShowerParticlePtr);

  /**
   *  Access to the original particle
   */
  inline PPtr original() const;

  /**
   *  Access to the colour isolated copy
   */
  inline PPtr copy() const;

  /**
   * Set the copy
   */
  inline void copy(PPtr);

  /**
   *  Whether the particle came from the hard process or was added by the matrix
   *  element correction
   */
  inline bool perturbative() const;

  /**
   *  Whether the particle came from the hard process or was added by the matrix
   *  element correction
   */
  inline void perturbative(bool);

  /**
   *  Access the \f$p_T\f$ of the hardest emission so far
   */
  inline Energy pT() const;

  /**
   *  Set the \f$p_T\f$ of the hardest emission so far
   */
  inline void pT(Energy);

  /**
   *  Has this particle radiated
   */
  inline bool hasEmitted() const;
  /**
   *  Set whether or not this particle has radiated
   */
  inline void hasEmitted(bool);

  /**
   *  The id of the particle
   */
  inline long id() const;

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
   *  Whether the particle came from the hard process or was added by the matrix
   *  element correction
   */
  bool _perturbative;

  /**
   *  Pointer to the ShowerParticle
   */
  ShowerParticlePtr _particle;

  /**
   *  \f$p_T\f$ of the current emmision
   */
  Energy _highestPt;

  /**
   *  Has there been radiation
   */
  bool _hasEmitted;

};
}

#include "ShowerProgenitor.icc"

#endif /* HERWIG_ShowerProgenitor_H */
