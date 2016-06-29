// -*- C++ -*-
//
// SplittingGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SplittingGenerator_H
#define HERWIG_SplittingGenerator_H
//
// This is the declaration of the SplittingGenerator class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/Shower/Base/Branching.h"
#include "Herwig/Shower/Base/SudakovFormFactor.h"
#include "SplittingGenerator.fh"
#include "Herwig/Shower/Base/ShowerParticle.h"
#include "Herwig/Shower/Base/ShowerKinematics.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 * 
 *  This class is responsible for creating, at the beginning of the Run,
 *  all the SplittingFunction objects and the corresponding
 *  SudakovFormFactor objects, and then of the generation of splittings 
 *  (radiation emissions) during the event.
 *  Many switches are defined in this class which allowed the user to turn on/off:
 *  -  each type of interaction (QCD, QED, EWK,...);
 *  - initial- and final-state radiation for all type of interactions;
 *  - initial- and final-state radiation for each type of interaction;
 *  - each type of splitting (\f$u\to ug\f$, \f$d\to dg\f$, \f$\ldots\f$, 
 *                            \f$g\to gg\f$, \f$g\to u\bar{u}\f$, \f$\ldots\f$).
 *
 *  These switches are useful mainly for debugging, but eventually can
 *  also be used for a "quick and dirty" estimation of systematic errors.
 *
 *  In the future it should be possible to implement in this class
 *
 *  -  the \f$1\to2\f$ azimuthal correlations for soft emission due to QCD coherence  
 *     using the ShowerParticle object provided in the input.
 *  -  Similarly having the \f$\rho-D\f$ matrix and the SplittingFunction pointer
 *     it should be possible to implement the spin correlations.
 *     
 *  @see SudakovFormFactor
 *  @see SplitFun
 *
 * @see \ref SplittingGeneratorInterfaces "The interfaces"
 * defined for SplittingGenerator.
 */
class SplittingGenerator: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SplittingGenerator() : _isr_Mode(1), _fsr_Mode(1), _deTuning(1.) {}
  //@}

public:

  /**
   *  Methods to select the next branching and reconstruct the kinematics
   */
  //@{
  /**
   * Choose a new forward branching for a time-like particle
   * The method returns:
   * - a pointer to a ShowerKinematics object, which 
   *     contains the information about the new scale and all other
   *     kinematics variables that need to be generated simultaneously;
   * - a pointer to the SudakovFormFactor object associated 
   *     with the chosen emission.
   * - The PDG codes of the particles in the branching,
   * as a Branching struct.
   *
   * In the case no branching has been generated, both the returned 
   * pointers are null ( ShoKinPtr() , tSudakovFFPtr() ).
   *
   * @param particle The particle to be evolved
   * @param enhance The factor by which to ehnace the emission of radiation
   * @param type The type of interaction to generate
   * @return The Branching struct for the branching
   */
  Branching chooseForwardBranching(ShowerParticle & particle,
				   double enhance,
				   ShowerInteraction::Type type) const; 

  /**
   * Select the next branching of a particles for the initial-state shower
   * in the particle's decay.
   * @param particle The particle being showerwed
   * @param maxscale The maximum scale
   * @param minmass Minimum mass of the particle after the branching
   * @param enhance The factor by which to ehnace the emission of radiation
   * @param type The type of interaction to generate
   * @return The Branching struct for the branching
   */
  Branching chooseDecayBranching(ShowerParticle & particle, 
				 const ShowerParticle::EvolutionScales & maxScales,
				 Energy minmass,double enhance,
				 ShowerInteraction::Type type) const; 

  /**
   * Choose a new backward branching for a space-like particle.
   * The method returns:
   * - a pointer to a ShowerKinematics object, which 
   *     contains the information about the new scale and all other
   *     kinematics variables that need to be generated simultaneously;
   * - a pointer to the SudakovFormFactor object associated 
   *     with the chosen emission.
   * - The PDG codes of the particles in the branching,
   * as a Branching struct.
   *
   * In the case no branching has been generated, both the returned 
   * pointers are null ( ShoKinPtr() , tSudakovFFPtr() ).
   *
   * @param particle The particle to be evolved
   * @param enhance The factor by which to ehnace the emission of radiation
   * @param beamparticle The beam particle
   * @param beam The BeamParticleData object
   * @param type The type of interaction to generate
   * @return The Branching struct for the branching
   */
  Branching 
  chooseBackwardBranching(ShowerParticle & particle,
			  PPtr beamparticle,
			  double enhance,
			  Ptr<BeamParticleData>::transient_const_pointer beam,
			  ShowerInteraction::Type type,
			  tcPDFPtr , Energy ) const;
  //@}

public:

  /**
   *  Access to the switches
   */
  //@{
  /**
   * It returns true/false if the initial-state radiation is on/off.
   */
  bool isISRadiationON() const { return _isr_Mode; }

  /**
   * It returns true/false if the final-state radiation is on/off.
   */
  bool isFSRadiationON() const { return _fsr_Mode; }
  //@}

  /**
   *  Methods to parse the information from the input files to create the 
   *  branchings
   */
  //@{
  /**
   *  Add a final-state splitting
   */
  string addFinalSplitting(string arg) { return addSplitting(arg,true); }

  /**
   *  Add an initial-state splitting
   */
  string addInitialSplitting(string arg) { return addSplitting(arg,false); }

  /**
   *  Add a final-state splitting
   */
  string deleteFinalSplitting(string arg) { return deleteSplitting(arg,true); }

  /**
   *  Add an initial-state splitting
   */
  string deleteInitialSplitting(string arg) { return deleteSplitting(arg,false); }
  //@}

  /**
   *  Access to the splittings
   */
  //@{
  /**
   *  Access the final-state branchings
   */
  const BranchingList & finalStateBranchings() const { return _fbranchings; }

  /**
   *  Access the initial-state branchings
   */
  const BranchingList & initialStateBranchings() const { return _bbranchings; }
  //@}

  /**
   * Set the factorization scale factor
   */
  void factorizationScaleFactor(double f);

  /**
   * Set the renormalization scale factor
   */
  void renormalizationScaleFactor(double f);

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans)
   ;

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

private:

  /**
   *  Add a branching to the map
   * @param ids PDG coeds of the particles in the branching
   * @param sudakov The SudakovFormFactor for the branching
   * @param final Whether this is an initial- or final-state branching 
   */
  void addToMap(const IdList & ids, const SudakovPtr & sudakov, bool final);

  /**
   * Remove a branching to the map
   * @param ids PDG coeds of the particles in the branching
   * @param sudakov The SudakovFormFactor for the branching
   * @param final Whether this is an initial- or final-state branching 
   */
  void deleteFromMap(const IdList & ids, const SudakovPtr & sudakov, bool final);

  /**
   * Obtain the reference vectors for a final-state particle
   * @param particle The particle
   * @param p The p reference vector
   * @param n The n reference vector
   */
  void finalStateBasisVectors(ShowerParticle particle, Lorentz5Momentum & p,
			      Lorentz5Momentum & n) const;

  /**
   * Add a splitting
   * @param in string to be parsed
   * @param final Whether this is an initial- or final-state branching 
   */
  string addSplitting(string in ,bool final);

  /**
   * Delete a splitting
   * @param in string to be parsed
   * @param final Whether this is an initial- or final-state branching 
   */
  string deleteSplitting(string in ,bool final);

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SplittingGenerator & operator=(const SplittingGenerator &);

private:

  /**
   *  Switches to control the radiation
   */
  //@{
  /**
   *  Is inqitial-state radiation on/off
   */
  bool _isr_Mode;

  /**
   *  Is final-state radiation on/off
   */
  bool _fsr_Mode;
  //@}

  /**
   *  List of the branchings and the appropriate Sudakovs for forward branchings
   */
  BranchingList _fbranchings;

  /**  
   * Lists of the branchings and the appropriate Sudakovs for backward branchings.
   */
  BranchingList _bbranchings;

  /**
   *   The detuning parameter
   */
  double _deTuning;
};

}

#endif /* HERWIG_SplittingGenerator_H */
