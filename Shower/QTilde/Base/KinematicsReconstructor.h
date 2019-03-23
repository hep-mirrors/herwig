// -*- C++ -*-
//
// KinematicsReconstructor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_KinematicsReconstructor_H
#define HERWIG_KinematicsReconstructor_H
//
// This is the declaration of the KinematicsReconstructor class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/Shower/Core/Base/ShowerParticle.h"
#include "Herwig/Shower/Core/Base/ShowerProgenitor.h"
#include "Herwig/Shower/Core/Base/ShowerTree.h"
#include "Herwig/Shower/Core/Base/HardTree.h"
#include "KinematicsReconstructor.fh"
#include <cassert>

namespace Herwig {

using namespace ThePEG;

  /**\ingroup Shower
   * Exception class
   * used to communicate failure of kinematics
   * reconstruction.
   */
  struct KinematicsReconstructionVeto {};


/** \ingroup Shower
 *
 * This class is responsible for the kinematical reconstruction 
 * after each showering step, and also for the necessary Lorentz boosts 
 * in order to preserve energy-momentum conservation in the overall collision,
 * and also the invariant mass and the rapidity of the hard subprocess system.
 * In the case of multi-step showering, there will be not unnecessary
 * kinematical reconstructions. 
 *
 * Notice:
 * - although we often use the term "jet" in either methods or variables names,
 *   or in comments, which could appear applicable only for QCD showering,
 *   there is indeed no "dynamics" represented in this class: only kinematics 
 *   is involved, as the name of this class remainds. Therefore it can be used
 *   for any kind of showers (QCD-,QED-,EWK-,... bremsstrahlung).
 * 
 * @see \ref KinematicsReconstructorInterfaces "The interfaces"
 * defined for KinematicsReconstructor.
 */
class KinematicsReconstructor: public Interfaced {

public:

  /**
   *  Methods to reconstruct the kinematics of a scattering or decay process
   */
  //@{
  /**
   * Given the ShowerTree for the shower from a hard process
   * the method does the reconstruction of the jets,
   * including the appropriate boosts (kinematics reshufflings)  
   * needed to conserve the total energy-momentum of the collision
   * and preserving the invariant mass and the rapidity of the 
   * hard subprocess system.
   */
  virtual bool reconstructHardJets(ShowerTreePtr hard,
				   const map<tShowerProgenitorPtr,
				   pair<Energy,double> > & pt,
				   ShowerInteraction type,
				   bool switchRecon) const=0;

  /**
   * Given the ShowerTree for a decay shower
   * the method does the reconstruction of the jets,
   * including the appropriate boosts (kinematics reshufflings)  
   * needed to conserve the total energy-momentum of the collision
   * and preserving the invariant mass and the rapidity of the 
   * hard subprocess system.
   */
  virtual bool reconstructDecayJets(ShowerTreePtr decay,
				    ShowerInteraction type) const=0;
  //@}

  /**
   *  Methods to invert the reconstruction of the shower for
   *  a scattering or decay process and calculate
   *  the variables used to generate the
   *  shower given the particles produced.
   *  This is needed for the CKKW and POWHEG approaches
   */
  //@{
  /**
   *  Given the particles, with a history which we wish to interpret
   *  as a shower reconstruct the variables used to generate the 
   * shower for a decay process
   */
  virtual bool deconstructDecayJets(HardTreePtr decay,ShowerInteraction) const=0;

  /**
   *  Given the particles, with a history which we wish to interpret
   *  as a shower reconstruct the variables used to generate the shower
   *  for a hard process
   */
  virtual bool deconstructHardJets(HardTreePtr hard,ShowerInteraction) const=0;
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KinematicsReconstructor & operator=(const KinematicsReconstructor &) = delete;
};

}

#endif /* HERWIG_KinematicsReconstructor_H */
