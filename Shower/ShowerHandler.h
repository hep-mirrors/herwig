// -*- C++ -*-
//
// ShowerHandler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerHandler_H
#define HERWIG_ShowerHandler_H
//
// This is the declaration of the ShowerHandler class.
//

#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ShowerVariation.h"
#include "Herwig/PDF/HwRemDecayer.fh"
#include "ThePEG/EventRecord/RemnantParticle.fh"
#include "UEBase.h"
#include "PerturbativeProcess.h"
#include "Herwig/MatrixElement/Matchbox/Matching/HardScaleProfile.h"
#include "ShowerHandler.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Shower
 *
 *  This class is the main driver of the shower: it is responsible for 
 *  the proper handling of all other specific collaborating classes
 *  and for the storing of the produced particles in the event record.
 * 
 *  @see \ref ShowerHandlerInterfaces "The interfaces"
 *
 *  @see ThePEG::CascadeHandler
 *  @see MPIHandler
 *  @see HwRemDecayer
 */
class ShowerHandler: public CascadeHandler {

public:

  /**
   * Typedef for a pair of ThePEG::RemnantParticle pointers.
   */
  typedef pair<tRemPPtr, tRemPPtr> RemPair;

public:

  /**
   *  Default constructor
   */
  ShowerHandler();

  /**
   *  Destructor
   */
  virtual ~ShowerHandler();

public:

  /**
   * The main method which manages the multiple interactions and starts
   * the shower by calling cascade(sub, lastXC).
   */
  virtual void cascade();

  /**
   *  pointer to "this", the current ShowerHandler.
   */
  static const tShowerHandlerPtr currentHandler() {
    assert(currentHandler_);
    return currentHandler_;
  }

public:

  /**
   * Hook to allow vetoing of event after showering hard sub-process
   * as in e.g. MLM merging.
   */
  virtual bool showerHardProcessVeto() const { return false; }

  /**
   * Return true, if this cascade handler will perform reshuffling from hard
   * process masses.
   */
  virtual bool isReshuffling() const { return true; }
  
  /**
   * Return true, if this cascade handler will put the final state
   * particles to their constituent mass. If false the nominal mass is used.
   */
  virtual bool retConstituentMasses() const { return useConstituentMasses_; }
  
  

  /**
   * Return true, if the shower handler can generate a truncated 
   * shower for POWHEG style events generated using Matchbox
   */
  virtual bool canHandleMatchboxTrunc() const { return false; }

  /**
   * Get the PDF freezing scale
   */
  Energy pdfFreezingScale() const { return pdfFreezingScale_; }
  
  /**
   * Get the local PDFs.
   */
  PDFPtr getPDFA() const {return PDFA_;}
  
  /**
   * Get the local PDFs.
   */
  PDFPtr getPDFB() const {return PDFB_;}
  
  /**
   * Return true if currently the primary subprocess is showered.
   */
  bool firstInteraction() const {
    if (!eventHandler()->currentCollision())return true;
    return ( subProcess_ == 
	     eventHandler()->currentCollision()->primarySubProcess() );
  }

  /**
   * Return the remnant decayer.
   */
  tHwRemDecPtr remnantDecayer() const { return remDec_; }

  /**
   *  Split the hard process into production and decays
   * @param tagged The tagged particles from the StepHandler
   * @param hard  The hard perturbative process
   * @param decay The decay particles
   */
  void splitHardProcess(tPVector tagged, PerturbativeProcessPtr & hard,
			DecayProcessMap & decay) const;

  /**
   * Information if the Showerhandler splits the hard process. 
   */
  bool doesSplitHardProcess()const {return splitHardProcess_;}

  /**
   *  Decay a particle.
   *  radPhotons switches the generation of photon
   *  radiation on/off.
   *  Required for Dipole Shower but not QTilde Shower.
   */
  tDMPtr decay(PerturbativeProcessPtr,
	       DecayProcessMap & decay,
	       bool radPhotons = false) const;

  
  /**
   * Cached lookup of decay modes.
   * Generator::findDecayMode() is not efficient.
   */
  tDMPtr findDecayMode(const string & tag) const;


  /**
   *  A struct to order the particles in the same way as in the DecayMode's
   */
  
  struct ParticleOrdering {

    bool operator() (tcPDPtr p1, tcPDPtr p2) const;

  };
  

  /**
   * A container for ordered particles required
   * for constructing tags for decay mode lookup.
   */
  typedef multiset<tcPDPtr,ParticleOrdering> OrderedParticles;


public:

  /**
   * @name Switches for initial- and final-state radiation
   */
  //@{
  /**
   *  Switch for any radiation
   */
  bool doRadiation() const {return doFSR_ || doISR_;}

  /**
   * Switch on or off final state radiation.
   */
  bool doFSR() const { return doFSR_;}
  
  /**
   * Switch on or off initial state radiation.
   */
  bool doISR() const { return doISR_;}
  //@}

public:

  /**
   * @name Switches for scales
   */
  //@{
  /**
   * Return true if maximum pt should be deduced from the factorization scale
   */
  bool hardScaleIsMuF() const { return maxPtIsMuF_; }

  /**
   * The factorization scale factor.
   */
  double factorizationScaleFactor() const {
    return factorizationScaleFactor_;
  }
  
  /**
   * The renormalization scale factor.
   */
  double renFac() const {
    return renormalizationScaleFactor_;
  }
  
  /**
   * The factorization scale factor.
   */
  double facFac() const {
    return factorizationScaleFactor_;
  }
  
  /**
   * The renormalization scale factor.
   */
  double renormalizationScaleFactor() const {
    return renormalizationScaleFactor_;
  }

  /**
   * The scale factor for the hard scale
   */
  double hardScaleFactor() const {
    return hardScaleFactor_;
  }
  /**
   * Return true, if the phase space restrictions of the dipole shower should
   * be applied.
   */
  bool restrictPhasespace() const { return restrictPhasespace_; }

  /**
   * Return profile scales
   */
  Ptr<HardScaleProfile>::tptr profileScales() const { return hardScaleProfile_; }

  /**
   * Return the relevant hard scale to be used in the profile scales
   */
  virtual Energy hardScale() const;

  /**
   * Return information about shower phase space choices
   */
  virtual int showerPhaseSpaceOption() const {
    assert(false && "not implemented in general");
    return -1;
  }
  //@}

public:

  /**
   * Access the shower variations
   */
  map<string,ShowerVariation>& showerVariations() {
    return showerVariations_;
  }

  /**
   * Return the shower variations
   */
  const map<string,ShowerVariation>& showerVariations() const {
    return showerVariations_;
  }

  /**
   * Access the current Weights
   */
  map<string,double>& currentWeights() {
    return currentWeights_;
  }

  /**
   * Return the current Weights
   */
  const map<string,double>& currentWeights() const {
    return currentWeights_;
  }

  /**
   * Change the current reweighting factor
   */
  void reweight(double w) {
    reweight_ = w;
  }

  /**
   * Return the current reweighting factor
   */
  double reweight() const {
    return reweight_;
  }

public:

  /** 
   * struct that is used to catch exceptions which are thrown
   * due to energy conservation issues of additional scatters
   */
  struct ExtraScatterVeto {};

  /** 
   * struct that is used to catch exceptions which are thrown
   * due to fact that the Shower has been invoked more than
   * a defined threshold on a certain configuration
   */
  struct ShowerTriesVeto {
    /** variable to store the number of attempts */
    const int tries;

    /** constructor */
    ShowerTriesVeto(int t) : tries(t) {}
  };

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

  /** @name Functions to perform the cascade
   */
  //@{
  /**
   * The main method which manages the showering of a subprocess.
   */
  virtual tPPair cascade(tSubProPtr sub, XCPtr xcomb);

  /**
   * Set up for the cascade
   */
  void prepareCascade(tSubProPtr sub) { 
    current_ = currentStep(); 
    subProcess_ = sub;
  } 

  /**
   *  Boost all the particles in the collision so that the collision always occurs
   * in the rest frame with the incoming particles along the z axis
   */
  void boostCollision(bool boost);
  //@}

protected:

  /**
   *  Set/unset the current shower handler
   */
  //@{
  /**
   *  Set the current handler
   */
  void setCurrentHandler() {
    currentHandler_ = tShowerHandlerPtr(this);
  }

  /**
   * Unset the current handler
   */
  void unSetCurrentHandler() {
    currentHandler_ = tShowerHandlerPtr();
  }
  //@}

protected:

  /**
   * @name Members relating to the underlying event and MPI
   */
  //@{
  /**
   * Return true if multiple parton interactions are switched on 
   * and can be used for this beam setup.
   */
  bool isMPIOn() const {
    return MPIHandler_ && MPIHandler_->beamOK();
  }

  /**
   * Access function for the MPIHandler, it should only be called after
   * checking with isMPIOn.
   */
  tUEBasePtr getMPIHandler() const {  
    assert(MPIHandler_);
    return MPIHandler_;    
  }

  /**
   *  Is a beam particle where hadronic structure is resolved
   */
  bool isResolvedHadron(tPPtr);

  /**
   * Get the remnants from the ThePEG::PartonBinInstance es and 
   * do some checks.
   */
  RemPair getRemnants(PBIPair incbins);

  /**
   *  Reset the PDF's after the hard collision has been showered
   */
  void setMPIPDFs();
  //@}

public:

  /**
   *  Check if a particle decays in the shower
   * @param id The PDG code for the particle
   */
  bool decaysInShower(long id) const {
    return ( particlesDecayInShower_.find( abs(id) ) != 
	     particlesDecayInShower_.end() );
  }

protected:

  /**
   *   Members to handle splitting up of hard process and decays
   */
  //@{
  /**
   *  Find decay products from the hard process and create decay processes
   * @param parent The parent particle
   * @param hard  The hard process
   * @param decay The decay processes
   */
  void findDecayProducts(PPtr parent, PerturbativeProcessPtr hard, DecayProcessMap & decay) const;

  /**
   * Find decay products from the hard process and create decay processes
   * @param parent The parent particle
   * @param hard  The parent hard process
   * @param decay The decay processes
   */
  void createDecayProcess(PPtr parent,PerturbativeProcessPtr hard, DecayProcessMap & decay) const;
  //@}

  /**
   * @name Functions to return information relevant to the process being showered
   */
  //@{
  /**
   * Return the currently used SubProcess.
   */
  tSubProPtr currentSubProcess() const {
    assert(subProcess_);
    return subProcess_;
  }

  /**
   *  Access to the incoming beam particles
   */
  tPPair incomingBeams() const {
    return incoming_;
  }
  //@}

protected:

  /**
   *  Weight handling for shower variations
   */
  //@
  /**
   * Combine the variation weights which have been encountered
   */
  void combineWeights();

 /**
   * Initialise the weights in currentEvent()
   */
  void initializeWeights();

  /**
   * Reset the current weights
   */
  void resetWeights();
  //@}

protected:

  /**
   * Return the maximum number of attempts for showering
   * a given subprocess.
   */
  unsigned int maxtry() const { return maxtry_; }

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

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerHandler & operator=(const ShowerHandler &) = delete;

private:

  /**
   *  pointer to "this", the current ShowerHandler.
   */
  static tShowerHandlerPtr  currentHandler_;

  /**
   * a MPIHandler to administer the creation of several (semihard) 
   * partonic interactions.
   */
  UEBasePtr MPIHandler_;

  /**
   *  Pointer to the HwRemDecayer
   */
  HwRemDecPtr remDec_;

private:

  /**
   *  Maximum tries for various stages of the showering process
   */
  //@{
  /**
   *  Maximum number of attempts for the
   *   main showering loop
   */
  unsigned int maxtry_;

  /**
   *  Maximum number of attempts for the regeneration of an additional
   *  scattering, before the number of scatters is reduced.
   */
  unsigned int maxtryMPI_;

  /**
   *  Maximum number of attempts for the regeneration of an additional
   *  hard scattering, before this event is vetoed.
   */
  unsigned int maxtryDP_;

  /**
   *  Maximum number of attempts to generate a decay
   */
  unsigned int maxtryDecay_;
  //@}

private:

  /**
   *  Factors for the various scales
   */
  //@{
  /**
   * The factorization scale factor.
   */
  double factorizationScaleFactor_;

  /**
   * The renormalization scale factor.
   */
  double renormalizationScaleFactor_;

  /**
   * The scale factor for the hard scale
   */
  double hardScaleFactor_;

  /**
   * True, if the phase space restrictions of the dipole shower should
   * be applied.
   */
  bool restrictPhasespace_;

  /**
   * True if maximum pt should be deduced from the factorization scale
   */
  bool maxPtIsMuF_;

  /**
   * The profile scales
   */
  Ptr<HardScaleProfile>::ptr hardScaleProfile_;
  //@}

private:

  /**
   *  Storage of information about the current event
   */
  //@{
  /**
   *  The incoming beam particles for the current collision
   */
  tPPair incoming_;

  /**
   *  Boost to get back to the lab
   */
  LorentzRotation boost_;

  /**
   *  Const pointer to the currently handeled ThePEG::SubProcess
   */
  tSubProPtr subProcess_;

  /**
   *  Const pointer to the current step
   */
  tcStepPtr current_;
  //@}

private:

  /**
   * PDFs to be used for the various stages and related parameters
   */
  //@{
  /**
   * The PDF freezing scale
   */
  Energy pdfFreezingScale_;

  /**
   * PDFs to be used for the various stages and related parameters
   */
  //@{
  /**
   * The PDF for beam particle A. Overrides the particle's own PDF setting.
   */
  PDFPtr PDFA_;

  /**
   * The PDF for beam particle B. Overrides the particle's own PDF setting.
   */
  PDFPtr PDFB_;

  /**
   * The PDF for beam particle A for remnant splitting. Overrides the particle's own PDF setting.
   */
  PDFPtr PDFARemnant_;

  /**
   * The PDF for beam particle B for remnant splitting. Overrides the particle's own PDF setting.
   */
  PDFPtr PDFBRemnant_;

  /**
   * The MPI PDF's to be used for secondary scatters.
   */
  pair <PDFPtr, PDFPtr> mpipdfs_;

  /**
   * The MPI PDF's to be used for secondary scatters.
   */
  pair <PDFPtr, PDFPtr> rempdfs_;

  /**
   * The MPI PDF's to be used for secondary scatters.
   */
  pair <PDFPtr, PDFPtr> remmpipdfs_;
  //@}

private:

  /**
   * @name Parameters for initial- and final-state radiation
   */
  //@{
  /**
   * Switch on or off final state radiation.
   */
  bool doFSR_;

  /**
   * Switch on or off initial state radiation.
   */
  bool doISR_;
  //@}

private:

  /**
   * @name Parameters for particle decays
   */
  //@{
  /**
   *  Whether or not to split into hard and decay trees
   */
  bool splitHardProcess_;

  /**
   *  PDG codes of the particles which decay during showering
   *  this is fast storage for use during running
   */
  set<long> particlesDecayInShower_;

  /**
   *  PDG codes of the particles which decay during showering
   *  this is a vector that is interfaced so they can be changed
   */
  vector<long> inputparticlesDecayInShower_;
  //@}

private:

  /**
   *  Parameters for the space-time model
   */
  //@{
  /**
   *   Whether or not to include spa-cetime distances in the shower
   */
  bool includeSpaceTime_;

  /**
   *  The minimum virtuality for the space-time model
   */
  Energy2 vMin_;
  //@}


private:
  
  /**
   *  Parameters for the constituent mass treatment.
   */
    //@{
  // True if shower should return constituent masses.
  bool useConstituentMasses_=true;
  //@}
private:

  /**
   *  Parameters relevant for reweight and variations
   */
  //@{
  /**
   * The shower variations
   */
  map<string,ShowerVariation> showerVariations_;

  /**
   * Command to add a shower variation
   */
  string doAddVariation(string);

  /**
   * A reweighting factor applied by the showering
   */
  double reweight_;

  /**
   * The shower variation weights
   */
  map<string,double> currentWeights_;
  //@}
};

}

#endif /* HERWIG_ShowerHandler_H */
