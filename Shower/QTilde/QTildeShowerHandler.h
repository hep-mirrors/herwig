// -*- C++ -*-
#ifndef Herwig_QTildeShowerHandler_H
#define Herwig_QTildeShowerHandler_H
//
// This is the declaration of the QTildeShowerHandler class.
//

#include "QTildeShowerHandler.fh"
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Shower/QTilde/SplittingFunctions/SplittingGenerator.h"
#include "Herwig/Shower/QTilde/Base/ShowerTree.h"
#include "Herwig/Shower/QTilde/Base/ShowerProgenitor.fh"
#include "Herwig/Shower/QTilde/Base/HardTree.h"
#include "Herwig/Shower/QTilde/Base/Branching.h"
#include "Herwig/Shower/QTilde/Base/ShowerVeto.h"
#include "Herwig/Shower/QTilde/Base/FullShowerVeto.h"
#include "Herwig/Shower/QTilde/Kinematics/KinematicsReconstructor.fh"
#include "Herwig/Shower/QTilde/Base/PartnerFinder.fh"
#include "Herwig/Shower/QTilde/SplittingFunctions/SudakovFormFactor.fh"
#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/Decay/HwDecayerBase.h"
#include "Herwig/MatrixElement/Matchbox/Matching/ShowerApproximation.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Utilities/Statistic.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The QTildeShowerHandler class.
 *
 * @see \ref QTildeShowerHandlerInterfaces "The interfaces"
 * defined for QTildeShowerHandler.
 */
class QTildeShowerHandler: public ShowerHandler {

public:
  
  /**
   *  Pointer to an XComb object
   */
  typedef Ptr<XComb>::pointer XCPtr;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  QTildeShowerHandler();

  /**
   * The destructor.
   */
  virtual ~QTildeShowerHandler();
  //@}

public:

  /**
   * At the end of the Showering, transform ShowerParticle objects
   * into ThePEG particles and fill the event record with them.
   * Notice that the parent/child relationships and the 
   * transformation from ShowerColourLine objects into ThePEG
   * ColourLine ones must be properly handled.
   */
  void fillEventRecord();
  
  /**
   * Return the relevant hard scale to be used in the profile scales
   */
  virtual Energy hardScale() const {
    return muPt;
  }

  /**
   * Hook to allow vetoing of event after showering hard sub-process
   * as in e.g. MLM merging.
   */
  virtual bool showerHardProcessVeto() const { return false; }
  
  /**
   *  Generate hard emissions for CKKW etc
   */
  virtual HardTreePtr generateCKKW(ShowerTreePtr tree) const;

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
   *  Get the SplittingGenerator
   */
  tSplittingGeneratorPtr splittingGenerator() const { return _splittingGenerator; }

  /**
   * Mode for hard emissions
   */
  int hardEmission() const {return _hardEmission;}
  //@}

  /**
   *  Connect the Hard and Shower trees
   */
  virtual void connectTrees(ShowerTreePtr showerTree, HardTreePtr hardTree, bool hard );

  /**
   *   Access to switches for spin correlations
   */
  //@{
  /**
   *  Soft correlations
   */
  unsigned int softCorrelations() const {
    return _softOpt;
  }

  /**
   *  Any correlations
   */
  virtual bool correlations() const {
    return spinCorrelations()!=0||_softOpt!=0;
  }
  //@}

public:
  /**
   *  Access methods to access the objects
   */
  //@{
  /**
   *  Access to the KinematicsReconstructor object
   */
  tKinematicsReconstructorPtr kinematicsReconstructor() const { return _reconstructor; }

  /**
   *  Access to the PartnerFinder object
   */
  tPartnerFinderPtr partnerFinder() const { return _partnerfinder; }

  //@}

protected:

  /**
   *   Perform the shower
   */
  void doShowering(bool hard,XCPtr);

  /**
   *  Generate the hard matrix element correction
   */
  virtual RealEmissionProcessPtr hardMatrixElementCorrection(bool);

  /**
   *  Generate the hardest emission
   */
  virtual void hardestEmission(bool hard);

  /**
   *  Set up for applying a matrix element correction
   */
  void setupMECorrection(RealEmissionProcessPtr real);

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
  virtual void setEvolutionPartners(bool hard,ShowerInteraction,
				    bool clear);

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
  virtual bool timeLikeShower(tShowerParticlePtr particle, ShowerInteraction,
			      Branching fb, bool first);

  /**
   * It does the backward evolution of the space-like input particle 
   * (and recursively for all its time-like radiation products).
   * accepting only emissions which conforms to the showerVariables.
   * If at least one emission has occurred then the method returns true
   * @param particle The particle to be showered
   * @param beam The beam particle
   */
  virtual bool spaceLikeShower(tShowerParticlePtr particle,PPtr beam,
			       ShowerInteraction); 

  /**
   * If does the forward evolution of the input on-shell particle
   * involved in a decay 
   * (and recursively for all its time-like radiation products).
   * accepting only emissions which conforms to the showerVariables.
   * @param particle    The particle to be showered
   * @param maxscale    The maximum scale for the shower.
   * @param minimumMass The minimum mass of the final-state system
   */
  virtual bool 
  spaceLikeDecayShower(tShowerParticlePtr particle,
		       const ShowerParticle::EvolutionScales & maxScales,
		       Energy minimumMass,ShowerInteraction,
		       Branching fb);

  /**
   * Truncated shower from a time-like particle
   */
  virtual bool truncatedTimeLikeShower(tShowerParticlePtr particle,
				       HardBranchingPtr branch,
				       ShowerInteraction type,
				       Branching fb, bool first);
 
  /**
   * Truncated shower from a space-like particle
   */
  virtual bool truncatedSpaceLikeShower(tShowerParticlePtr particle,PPtr beam,
					HardBranchingPtr branch,
					ShowerInteraction type);

  /**
   * Truncated shower from a time-like particle
   */
  virtual bool truncatedSpaceLikeDecayShower(tShowerParticlePtr particle,
					     const ShowerParticle::EvolutionScales & maxScales,
					     Energy minimumMass, HardBranchingPtr branch,
					     ShowerInteraction type, Branching fb);
  //@}

  /**
   *  Switches for matrix element corrections
   */
  //@{
  /**
   * Any ME correction?   
   */
  bool MECOn() const { 
    return _hardEmission == 1;
  }

  /**
   * Any hard ME correction? 
   */
  bool hardMEC() const {
    return _hardEmission == 1 && (_meCorrMode == 1 || _meCorrMode == 2);
  }

  /**
   * Any soft ME correction? 
   */
  bool softMEC() const {
    return _hardEmission == 1 && (_meCorrMode == 1 || _meCorrMode > 2);
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
   * Set/Get the current tree being evolver for inheriting classes
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
  void setupMaximumScales(const vector<ShowerProgenitorPtr> &,XCPtr);

  /**
   *  find the relevant hard scales for profile scales. 
   */
  void setupHardScales(const vector<ShowerProgenitorPtr> &,XCPtr);

  /**
   *  Convert the HardTree into an extra shower emission 
   */
  void convertHardTree(bool hard,ShowerInteraction type);

protected:

  /**
   * Find the parton extracted from the incoming particle after ISR
   */
  PPtr findFirstParton(tPPtr seed) const;

  /**
   * Fix Remnant connections after ISR
   */
  tPPair remakeRemnant(tPPair oldp); 

protected:

  /**
   *  Start the shower of a timelike particle
   */
  virtual bool startTimeLikeShower(ShowerInteraction);

  /**
   *  Update of the time-like stuff
   */
  void updateHistory(tShowerParticlePtr particle);

  /**
   *  Start the shower of a spacelike particle
   */
  virtual bool startSpaceLikeShower(PPtr,ShowerInteraction);

  /**
   *  Start the shower of a spacelike particle
   */
  virtual bool 
  startSpaceLikeDecayShower(const ShowerParticle::EvolutionScales & maxScales,
			    Energy minimumMass,ShowerInteraction);

  /**
   * Select the branching for the next time-like emission
   */
  Branching selectTimeLikeBranching(tShowerParticlePtr particle,
				    ShowerInteraction type,
				    HardBranchingPtr branch);

  /**
   * Select the branching for the next space-like emission in a decay
   */
  Branching selectSpaceLikeDecayBranching(tShowerParticlePtr particle,
					  const ShowerParticle::EvolutionScales & maxScales,
					  Energy minmass,ShowerInteraction type,
					  HardBranchingPtr branch);
  /**
   *  Create the timelike child of a branching
   */
  ShowerParticleVector createTimeLikeChildren(tShowerParticlePtr particle,
					      IdList ids);

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
  bool hardOnly() const {return _limitEmissions==3;}

  /**
   *  Check the flags
   */
  void checkFlags();

  /**
   *
   */
  void addFSRUsingDecayPOWHEG(HardTreePtr ISRTree);

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
   * The main method which manages the showering of a subprocess.
   */
  virtual tPPair cascade(tSubProPtr sub, XCPtr xcomb);

  /**
   *  Decay a ShowerTree
   */
  void decay(ShowerTreePtr tree, ShowerDecayMap & decay);

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

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeShowerHandler & operator=(const QTildeShowerHandler &) = delete;

private:


  /**
   *  Stuff from the ShowerHandler
   */
  //@{

  /**
   *  The ShowerTree for the hard process
   */
  ShowerTreePtr hard_;

  /**
   *  The ShowerTree for the decays
   */
  ShowerDecayMap decay_;

  /**
   *  The ShowerTrees for which the initial shower 
   */
  vector<ShowerTreePtr> done_;
  //@}

private :

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
   *  Control of the reconstruction option
   */
  unsigned int _evolutionScheme;

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
   *  Full Shower Vetoes
   */
  vector<FullShowerVetoPtr> _fullShowerVetoes;

  /**
   *  Number of iterations for reweighting
   */
  unsigned int _nReWeight;

  /**
   *  Whether or not we are reweighting
   */
  bool _reWeight;

  /**
   *  number of IS emissions
   */
  unsigned int _nis;

  /**
   *  Number of FS emissions
   */
  unsigned int _nfs;

  /**
   *  The option for wqhich interactions to use
   */
  ShowerInteraction interaction_;

  /**
   *  Truncated shower switch
   */
  bool _trunc_Mode;

  /**
   *  Mode for the hard emissions
   */
  int _hardEmission;

  /**
   *  Option for the kernal for soft correlations
   */
  unsigned int _softOpt;

  /**
   *  Option for hard radiation in POWHEG events
   */
  bool _hardPOWHEG;

  /**
   * True if no warnings about incorrect hard emission
   * mode setting have been issued yet
   */
  static bool _hardEmissionWarn;

  /**
   * True if no warnings about missing truncated shower 
   * have been issued yet
   */
  static bool _missingTruncWarn;

  /**
   * The relevant hard scale to be used in the profile scales
   */
  Energy muPt;

private:
  /**
   *  Pointer to the various objects
   */
  //@{
  /**
   *  Pointer to the KinematicsReconstructor object
   */
  KinematicsReconstructorPtr _reconstructor;

  /**
   *  Pointer to the PartnerFinder object
   */
  PartnerFinderPtr _partnerfinder;

  //@}

};

}

#endif /* HERWIG_QTildeShowerHandler_H */
