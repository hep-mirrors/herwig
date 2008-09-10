// -*- C++ -*-
//
// SplittingGenerator.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SplittingGenerator_H
#define HERWIG_SplittingGenerator_H
//
// This is the declaration of the SplittingGenerator class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower/Base/SudakovFormFactor.h"
#include "SplittingGenerator.fh"
#include "Herwig++/Shower/Base/ShowerKinematics.h"
#include "ThePEG/Utilities/Rebinder.h"
#include<vector>     

namespace Herwig {

using namespace ThePEG;

/**
 *  Forward declaration of the ShowerParticle class
 */
class ShowerParticle;

/**
 *  typedef to pair the SudakovFormFactor and the particles in a branching
 */
typedef pair<SudakovPtr,IdList> BranchingElement;

/**
 *  typedef to pair the PDG code of the particle and the BranchingElement
 */
typedef multimap<long,BranchingElement> BranchingList;

/**
 *  typedef to create a structure which can be inserted into a BranchingList
 */
typedef pair<long, BranchingElement> BranchingInsert; 

/** \ingroup Shower
 *  The branching struct is used to store information on the branching.
 *  The kinematics variable is a pointer to the ShowerKinematics for the branching
 *  The sudakov variable is a pointer to the SudakovFormFactor for the branching
 *  The ids  variable is the list of particles in the branching
 */
struct Branching {
 
  /**
   *  Pointer to the ShowerKinematics object for the branching
   */
  ShoKinPtr kinematics;
  
  /**
   *  PDG codes of the particles in the branching
   */
  IdList ids; 

  /**
   *  The SudakovFormFactor for the branching
   */
  tSudakovPtr sudakov;

  /**
   *  Constructor for the struct
   * @param a pointer to the ShowerKinematics object for the branching
   * @param c PDG codes of the particles in the branching
   * @param d The SudakovFormFactor for the branching
   */
  Branching(ShoKinPtr a, IdList c,tSudakovPtr d) : kinematics(a), ids(c), sudakov(d) {}

  /**
   *  Default constructor
   */
  Branching() {}
};

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
  SplittingGenerator() : _isr_Mode(1), _fsr_Mode(1) {}
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
   * @return The Branching struct for the branching
   */
  Branching chooseForwardBranching(ShowerParticle & particle,
				   double enhance) const; 

  /**
   * Select the next branching of a particles for the initial-state shower
   * in the particle's decay.
   * @param particle The particle being showerwed
   * @param maxscale The maximum scale
   * @param minmass Minimum mass of the particle after the branching
   * @param enhance The factor by which to ehnace the emission of radiation
   * @return The Branching struct for the branching
   */
  Branching chooseDecayBranching(ShowerParticle & particle, 
				 Energy maxscale,
				 Energy minmass,double enhance) const; 

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
   * @param beam The beam particle
   * @return The Branching struct for the branching
   */
  Branching 
  chooseBackwardBranching(ShowerParticle & particle,
			  PPtr beam,
			  double enhance,
			  Ptr<BeamParticleData>::transient_const_pointer) const;
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
    throw(RebindException);

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
   *  Obtain the reference vectors for a final-state particle
   * @param particle The particle
   * @param type The  type of the interaction
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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SplittingGenerator> initSplittingGenerator;

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
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SplittingGenerator. */
template <>
struct BaseClassTrait<Herwig::SplittingGenerator,1> {
  /** Typedef of the first base class of SplittingGenerator. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SplittingGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SplittingGenerator>
  : public ClassTraitsBase<Herwig::SplittingGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SplittingGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SplittingGenerator is implemented. It may also include several, space-separated,
   * libraries if the class SplittingGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SplittingGenerator_H */
