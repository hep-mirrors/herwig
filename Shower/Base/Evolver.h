// -*- C++ -*-
//
// Evolver.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Evolver_H
#define HERWIG_Evolver_H
//
// This is the declaration of the Evolver class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingGenerator.h"
#include "ShowerModel.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ShowerTree.h"
#include "ShowerProgenitor.fh"
#include "Herwig++/Shower/ShowerHandler.fh"
#include "Branching.h"
#include "ShowerVeto.h"
#include "HardTree.h"
#include "ThePEG/Handlers/XComb.h"
#include "Evolver.fh"
#include "Herwig++/MatrixElement/HwMEBase.h"
#include "Herwig++/Decay/HwDecayerBase.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 * The Evolver class class performs the sohwer evolution of hard scattering 
 * and decay processes in Herwig++.
 *
 * @see \ref EvolverInterfaces "The interfaces"
 * defined for Evolver.
 */
class Evolver: public Interfaced {

/**
 *  The ShowerHandler is a friend to set some parameters at initialisation
 */
friend class ShowerHandler;

public:
  
  /**
   *  Pointer to an XComb object
   */
  typedef Ptr<XComb>::pointer XCPtr;

public:

  /**
   *  Default Constructor
   */
  Evolver() : _maxtry(100), _meCorrMode(1), _hardVetoMode(1),
	      _hardVetoRead(0), _reconOpt(0), _hardVetoReadOption(false),
	      _iptrms(ZERO), _beta(0.), _gamma(ZERO), _iptmax(),
	      _limitEmissions(0), _initialenhance(1.), _finalenhance(1.),
	      _hardonly(false), _trunc_Mode(true), _hardEmissionMode(0),
	      _colourEvolutionMethod(0)
  {}

  /**
   *  Members to perform the shower
   */
  //@{
  /**
   * Perform the shower of the hard process
   */
  virtual void showerHardProcess(ShowerTreePtr,XCPtr);

  /**
   * Perform the shower of a decay
   */
  virtual void showerDecay(ShowerTreePtr);
  //@}

  /**
   *  Access to the flags and shower variables
   */
  //@{
  /**
   *  Is there any showering switched on
   */
  bool showeringON() const { return isISRadiationON() || isFSRadiationON(); }

  /**
   * It returns true/false if the initial-state radiation is on/off.
   */
  bool isISRadiationON() const { return _splittingGenerator->isISRadiationON(); }

  /**
   * It returns true/false if the final-state radiation is on/off.
   */
  bool isFSRadiationON() const { return _splittingGenerator->isFSRadiationON(); }

  /**
   *  Get the ShowerModel
   */ 
  ShowerModelPtr showerModel() const {return _model;}

  /**
   *  Get the SplittingGenerator
   */
  tSplittingGeneratorPtr splittingGenerator() const { return _splittingGenerator; }
  //@}

  /**
   *  Connect the Hard and Shower trees
   */
  virtual void connectTrees(ShowerTreePtr showerTree, HardTreePtr hardTree, bool hard )const;

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

  /**
   *  Generate the hard matrix element correction
   */
  virtual void hardMatrixElementCorrection(bool);

  /**
   *  Generate the hardest emission
   */
  virtual void hardestEmission(bool hard);

  /**
   * Extract the particles to be showered, set the evolution scales
   * and apply the hard matrix element correction
   * @param hard Whether this is a hard process or decay
   * @return The particles to be showered
   */
  virtual vector<ShowerProgenitorPtr> setupShower(bool hard);

  /**
   *  set the colour partners
   */
  virtual void setEvolutionPartners(bool hard,ShowerInteraction::Type);

  /**
   *  Methods to perform the evolution of an individual particle, including
   *  recursive calling on the products
   */
  //@{
  /**
   * It does the forward evolution of the time-like input particle
   * (and recursively for all its radiation products).
   * accepting only emissions which conforms to the showerVariables
   * and soft matrix element correction.
   * If at least one emission has occurred then the method returns true.
   * @param particle The particle to be showered
   */
  virtual bool timeLikeShower(tShowerParticlePtr particle, ShowerInteraction::Type,
			      bool first); 

  /**
   * It does the backward evolution of the space-like input particle 
   * (and recursively for all its time-like radiation products).
   * accepting only emissions which conforms to the showerVariables.
   * If at least one emission has occurred then the method returns true
   * @param particle The particle to be showered
   * @param beam The beam particle
   */
  virtual bool spaceLikeShower(tShowerParticlePtr particle,PPtr beam,
			       ShowerInteraction::Type); 

  /**
   * If does the forward evolution of the input on-shell particle
   * involved in a decay 
   * (and recursively for all its time-like radiation products).
   * accepting only emissions which conforms to the showerVariables.
   * @param particle    The particle to be showered
   * @param maxscale    The maximum scale for the shower.
   * @param minimumMass The minimum mass of the final-state system
   */
  virtual bool spaceLikeDecayShower(tShowerParticlePtr particle,
				    Energy maxscale,
				    Energy minimumMass,
				    ShowerInteraction::Type);

  /**
   * Truncated shower from a time-like particle
   */
  virtual bool truncatedTimeLikeShower(tShowerParticlePtr particle,
				       HardBranchingPtr branch,
				       ShowerInteraction::Type type);
 
  /**
   * Truncated shower from a space-like particle
   */
  virtual bool truncatedSpaceLikeShower(tShowerParticlePtr particle,PPtr beam,
					HardBranchingPtr branch,
					ShowerInteraction::Type type);

  /**
   * Truncated shower from a time-like particle
   */
  virtual bool truncatedSpaceLikeDecayShower(tShowerParticlePtr particle,
					     Energy maxscale, Energy minimumMass,
					     HardBranchingPtr branch,
					     ShowerInteraction::Type type);
  //@}

  /**
   *  Switches for matrix element corrections
   */
  //@{
  /**
   * Any ME correction?   
   */
  bool MECOn() const {
    return _meCorrMode > 0 && _hardEmissionMode==0;
  }

  /**
   * Any hard ME correction? 
   */
  bool hardMEC() const {
    return (_meCorrMode == 1 || _meCorrMode == 2) && _hardEmissionMode==0;
  }

  /**
   * Any soft ME correction? 
   */
  bool softMEC() const {
    return (_meCorrMode == 1 || _meCorrMode > 2) && _hardEmissionMode==0;
  }
  //@}

  /**
   * Is the truncated shower on?
   */
  bool isTruncatedShowerON() const {return _trunc_Mode;}

  /**
   *  Switch for intrinsic pT
   */
  //@{
  /**
   * Any intrinsic pT?
   */
  bool ipTon() const {
    return _iptrms != ZERO || ( _beta == 1.0 && _gamma != ZERO && _iptmax !=ZERO );
  }
   //@}  

  /**@name Additional shower vetoes */
  //@{
  /**
   * Insert a veto.
   */
  void addVeto (ShowerVetoPtr v) { _vetoes.push_back(v); }

  /**
   * Remove a veto.
   */
  void removeVeto (ShowerVetoPtr v) { 
    vector<ShowerVetoPtr>::iterator vit = find(_vetoes.begin(),_vetoes.end(),v);
    if (vit != _vetoes.end())
      _vetoes.erase(vit);
  }

  //@}

  /**
   *  Switches for vetoing hard emissions
   */
  //@{
  /**
   * Vetos on? 
   */
  bool hardVetoOn() const { return _hardVetoMode > 0; }

  /**
   * veto hard emissions in IS shower?
   */
  bool hardVetoIS() const { return _hardVetoMode == 1 || _hardVetoMode == 2; }

  /**
   * veto hard emissions in FS shower?
   */
  bool hardVetoFS() const { return _hardVetoMode == 1 || _hardVetoMode > 2; }

  /**
   * veto hard emissions according to lastScale from XComb? 
   */
  bool hardVetoXComb() const {return (_hardVetoRead == 1);}

  /**
   * Returns true if the hard veto read-in is to be applied to only
   * the primary collision and false otherwise.
   */
  bool hardVetoReadOption() const {return _hardVetoReadOption;}
  //@}

  /**
   *  Enhancement factors for radiation needed to generate the soft matrix
   *  element correction.
   */
  //@{
  /**
   *  Access the enhancement factor for initial-state radiation
   */
  double initialStateRadiationEnhancementFactor() const { return _initialenhance; }

  /**
   *  Access the enhancement factor for final-state radiation
   */
  double finalStateRadiationEnhancementFactor() const { return _finalenhance; }

  /**
   *  Set the enhancement factor for initial-state radiation
   */
  void initialStateRadiationEnhancementFactor(double in) { _initialenhance=in; }

  /**
   *  Set the enhancement factor for final-state radiation
   */
  void finalStateRadiationEnhancementFactor(double in) { _finalenhance=in; }
  //@}

  /**
   *  Access to set/get the HardTree currently beinging showered
   */
  //@{
  /**
   *  The HardTree currently being showered
   */
  tHardTreePtr hardTree() {return _hardtree;}

  /**
   *  The HardTree currently being showered
   */
  void hardTree(tHardTreePtr in) {_hardtree = in;}
  //@}

  /**
   * Access/set the beam particle for the current initial-state shower
   */
  //@{
  /**
   *  Get the beam particle data
   */
  Ptr<BeamParticleData>::const_pointer beamParticle() const { return _beam; }

  /**
   *  Set the beam particle data
   */
  void setBeamParticle(Ptr<BeamParticleData>::const_pointer in) { _beam=in; }
  //@}

  /**
   * Set/Get the current tree being evolverd for inheriting classes
   */
  //@{
  /**
   * Get the tree
   */
  tShowerTreePtr currentTree() { return _currenttree; }

  /**
   * Set the tree
   */
  void currentTree(tShowerTreePtr tree) { _currenttree=tree; }

  //@}

  /**
   *  Access the maximum number of attempts to generate the shower
   */
  unsigned int maximumTries() const { return _maxtry; }

  /**
   * Set/Get the ShowerProgenitor for the current shower
   */
  //@{
  /**
   *  Access the progenitor
   */
  ShowerProgenitorPtr progenitor() { return _progenitor; }

  /**
   *  Set the progenitor
   */
  void progenitor(ShowerProgenitorPtr in) { _progenitor=in; }
  //@}

  /**
   *  Calculate the intrinsic \f$p_T\f$.
   */
  virtual void generateIntrinsicpT(vector<ShowerProgenitorPtr>);

  /**
   *  Access to the intrinsic \f$p_T\f$ for inheriting classes
   */
  map<tShowerProgenitorPtr,pair<Energy,double> > & intrinsicpT() { return _intrinsic; }

  /**
   *  find the maximally allowed pt acc to the hard process. 
   */
  void setupMaximumScales(ShowerTreePtr, const vector<ShowerProgenitorPtr> &,XCPtr);

protected:

  /**
   *  Start the shower of a timelike particle
   */
  virtual bool startTimeLikeShower(ShowerInteraction::Type);

  /**
   *  Update of the time-like stuff
   */
  void updateHistory(tShowerParticlePtr particle);

  /**
   *  Start the shower of a spacelike particle
   */
  virtual bool startSpaceLikeShower(PPtr,ShowerInteraction::Type);

  /**
   *  Start the shower of a spacelike particle
   */
  virtual bool startSpaceLikeDecayShower(Energy maxscale,Energy minimumMass,
					 ShowerInteraction::Type);

  /**
   *  Vetos for the timelike shower
   */
  virtual bool timeLikeVetoed(const Branching &,ShowerParticlePtr);

  /**
   *  Vetos for the spacelike shower
   */
  virtual bool spaceLikeVetoed(const Branching &,ShowerParticlePtr);

  /**
   *  Vetos for the spacelike shower
   */
  virtual bool spaceLikeDecayVetoed(const Branching &,ShowerParticlePtr);

  /**
   *  Only generate the hard emission, for testing only.
   */
  bool hardOnly() const {return _hardonly;}

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

private:

  /**
   * Get the octet -> octet octet reduction factor.
   */
  double getReductionFactor(tShowerParticlePtr particle);

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<Evolver> initEvolver;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Evolver & operator=(const Evolver &);

private:

  /**
   *  Pointer to the model for the shower evolution model
   */
  ShowerModelPtr _model;

  /**
   * Pointer to the splitting generator
   */
  SplittingGeneratorPtr _splittingGenerator;

  /**
   *  Maximum number of tries to generate the shower of a particular tree
   */
  unsigned int _maxtry;

  /**
   * Matrix element correction switch
   */
  unsigned int _meCorrMode; 

  /**
   * Hard emission veto switch
   */
  unsigned int _hardVetoMode; 

  /**
   * Hard veto to be read switch
   */
  unsigned int _hardVetoRead; 

  /**
   *  Control of the reconstruction option
   */
  unsigned int _reconOpt;

  /**
   * If hard veto pT scale is being read-in this determines
   * whether the read-in value is applied to primary and 
   * secondary (MPI) scatters or just the primary one, with
   * the usual computation of the veto being performed for
   * the secondary (MPI) scatters.
   */
  bool _hardVetoReadOption; 

  /**
   * rms intrinsic pT of Gaussian distribution
   */
  Energy _iptrms;

   /**
   * Proportion of inverse quadratic intrinsic pT distribution
   */
  double _beta;

  /**
   * Parameter for inverse quadratic: 2*Beta*Gamma/(sqr(Gamma)+sqr(intrinsicpT))
   */
  Energy _gamma;

  /**
   * Upper bound on intrinsic pT for inverse quadratic
   */
  Energy _iptmax;

  /**
   *  Limit the number of emissions for testing
   */
  unsigned int _limitEmissions;
  
  /**
   *  The progenitor of the current shower
   */
  ShowerProgenitorPtr _progenitor;

  /**
   *  Matrix element
   */
  HwMEBasePtr _hardme;

  /**
   *  Decayer
   */
  HwDecayerBasePtr _decayme;

  /**
   * The ShowerTree currently being showered
   */
  ShowerTreePtr _currenttree;

  /**
   *  The HardTree currently being showered
   */
  HardTreePtr _hardtree;

  /**
   *  Radiation enhancement factors for use with the veto algorithm
   *  if needed by the soft matrix element correction 
   */
  //@{
  /**
   *  Enhancement factor for initial-state radiation
   */
  double _initialenhance;

  /**
   *  Enhancement factor for final-state radiation
   */
  double _finalenhance;
  //@}

  /**
   *  The beam particle data for the current initial-state shower
   */
  Ptr<BeamParticleData>::const_pointer _beam;

  /**
   *  Storage of the intrinsic \f$p_t\f$ of the particles
   */
  map<tShowerProgenitorPtr,pair<Energy,double> > _intrinsic;

  /**
   * Vetoes
   */
  vector<ShowerVetoPtr> _vetoes;

  /**
   *  number of IS emissions
   */
  unsigned int _nis;

  /**
   *  Number of FS emissions
   */
  unsigned int _nfs;

  /**
   *  Only generate the emission from the hardest emission
   *  generate for testing only
   */
  bool _hardonly;

 /**
   *  Truncated shower switch
   */
  bool _trunc_Mode;
  
  /**
   *  Count of the number of truncated emissions
   */
  unsigned int _truncEmissions;

  /**
   *  Mode for the hard emissions
   */
  unsigned int _hardEmissionMode;

  /**
   * Colour evolution method
   */
  int _colourEvolutionMethod;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Evolver. */
template <>
struct BaseClassTrait<Herwig::Evolver,1> {
  /** Typedef of the first base class of Evolver. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Evolver class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Evolver>
  : public ClassTraitsBase<Herwig::Evolver> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::Evolver"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Evolver is implemented. It may also include several, space-separated,
   * libraries if the class Evolver depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_Evolver_H */
