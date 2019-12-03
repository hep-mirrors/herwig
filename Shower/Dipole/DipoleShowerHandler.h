// -*- C++ -*-
//
// DipoleShowerHandler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleShowerHandler_H
#define HERWIG_DipoleShowerHandler_H
//
// This is the declaration of the DipoleShowerHandler class.
//

#include "Herwig/Shower/ShowerHandler.h"

#include "Herwig/Shower/Dipole/DipoleShowerHandler.fh"
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"
#include "Herwig/Shower/Dipole/Base/DipoleSplittingReweight.h"
#include "Herwig/Shower/Dipole/Kernels/DipoleSplittingKernel.h"
#include "Herwig/Shower/Dipole/Base/DipoleSplittingGenerator.h"
#include "Herwig/Shower/Dipole/Base/DipoleEventRecord.h"
#include "Herwig/Shower/Dipole/Base/DipoleEvolutionOrdering.h"
#include "Herwig/Shower/Dipole/Base/DipoleEventReweight.h"
#include "Herwig/Shower/Dipole/Utility/ConstituentReshuffler.h"
#include "Herwig/Shower/Dipole/Utility/IntrinsicPtGenerator.h"
#include "Herwig/MatrixElement/Matchbox/Base/MergerBase.h"
#include "Herwig/MatrixElement/Matchbox/Matching/ShowerApproximation.h"

#include "Herwig/Shower/Dipole/SpinCorrelations/DipoleVertexRecord.h"
#include "Herwig/MatrixElement/Matchbox/Utility/DensityOperator.h"

#include <tuple>

namespace Herwig {

using namespace ThePEG;

/** 
 * \ingroup DipoleShower
 * \author Simon Platzer, Stephen Webster
 *
 * \brief The DipoleShowerHandler class manages the showering using
 * the dipole shower algorithm.
 *
 * @see \ref DipoleShowerHandlerInterfaces "The interfaces"
 * defined for DipoleShowerHandler.
 */
class DipoleShowerHandler: public ShowerHandler {


 friend class Merger;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleShowerHandler();

  /**
   * The destructor.
   */
  virtual ~DipoleShowerHandler();
  //@}

public:


  inline void colourPrint();

  /**
   * Indicate a problem in the shower.
   */
  struct RedoShower {};

  /**
   * Insert an additional splitting kernel.
   */
  void addSplitting(Ptr<DipoleSplittingKernel>::ptr sp) {
    kernels.push_back(sp);
  }

  /**
   * Reset the alpha_s for all splitting kernels.
   */
  void resetAlphaS(Ptr<AlphaSBase>::tptr);
  
  
  virtual void cascade(tPVector); 

  /**
   * Reset the splitting reweight for all splitting kernels.
   */
  void resetReweight(Ptr<DipoleSplittingReweight>::tptr);

  /**
   * Return true, if the shower handler can generate a truncated 
   * shower for POWHEG style events generated using Matchbox
   */
  virtual bool canHandleMatchboxTrunc() const { return false; }

  /**
   * Return true, if this cascade handler will perform reshuffling from hard
   * process masses.
   */
  virtual bool isReshuffling() const { return false; }

  /**
   * Return the relevant hard scale to be used in the profile scales
   */
  virtual Energy hardScale() const {
    return muPt;
  }
  
  /**
   * Calculate the alpha_s value the shower uses for Q.
   */
  
  double as(Energy Q)const{return theGlobalAlphaS->value(sqr(Q));}
  
  /**
   * Return the number of scale dependent active flavours from 
   * the alpha_s object.
   */

  double Nf(Energy Q)const{return theGlobalAlphaS->Nf(sqr(Q));}

  /*
   * Access the vertex record
   */
  DipoleVertexRecord& vertexRecord() { return theVertexRecord; }

  /**
   * Return the event record
   */
  const DipoleVertexRecord& vertexRecord() const { return theVertexRecord; }


  /**
   * Set the pointer to the Merging Helper.
   * Used by the merging factory.
   */
  void setMerger(Ptr<MergerBase>::ptr mh){theMergingHelper=mh;}



public:

  /**
   * Return the dictionary for the incoming particles and outgoing partons 
   * in the event, will be empty if subleading Nc is turned off.
   */
  const map<PPtr,size_t>& particleIndices() const { return eventRecord().particleIndices(); }

  /**
   * Return the particle data vector of the particles in the event, will 
   * be empty if subleading Nc is turned off.
   */
  const cPDVector& particlesAfter() const { return eventRecord().particlesAfter(); }

  /**
   * Return the density operator from the event record.
   */
  DensityOperator& densityOperator() { return eventRecord().densityOperator(); }

  /**
   * Return the colour matrix element correction map, will be empty if
   * subleading Nc is turned off. NOT USED, REMOVE
   */
  const map<pair<vector<PDT::Colour>,pair<size_t,size_t> >,double>& correlatorMap() const {
    return eventRecord().densityOperator().correlatorMap();
  }

  /**
   * Return the colour basis used in the density operator. NOT USED, REMOVE
   */
  Ptr<ColourBasis>::tptr colourBasis() {
    return eventRecord().densityOperator().colourBasis();
  }

  /**
   * Return the continue subleading Nc flag from the event record.
   */
  bool continueSubleadingNc() const { return eventRecord().getContinueSubleadingNc(); }

protected:

  /**
   * Add the possible splitting candidates given a pair of emitting particles
   */
  void addCandidates(PPair particles, list<DipoleSplittingInfo>& clist) const;

  /**
   * Get all possible radiating dipoles
   */
  void getCandidates(list<DipoleSplittingInfo>& clist) const;

  /**
   * Perform a splitting independent of any chains
   */
  void performSplitting(DipoleSplittingInfo&) const;

  /**
   * Generate the next subleading Nc improved splitting and return the scale
   */
  Energy nextSubleadingSplitting(Energy hardPt,
				 Energy optHardPt, Energy optCutoff,
				 const bool decay);

protected:

  typedef multimap<DipoleIndex,Ptr<DipoleSplittingGenerator>::ptr> GeneratorMap;

  /**
   * The main method which manages the showering of a subprocess.
   */
  virtual tPPair cascade(tSubProPtr sub, XCombPtr xcomb) {
    return cascade(sub,xcomb,ZERO,ZERO);
  }

  /**
   * The main method which manages the showering of a subprocess.
   */
  tPPair cascade(tSubProPtr sub, XCombPtr xcomb, 
		 Energy optHardPt, Energy optCutoff);

  /**
   * Build splitting generators for the given
   * dipole index.
   */
  void getGenerators(const DipoleIndex&,
		     Ptr<DipoleSplittingReweight>::tptr rw =
  		     Ptr<DipoleSplittingReweight>::tptr());

  /**
   * Setup the hard scales.
   */
  void hardScales(Energy2 scale);

  /**
   * Setup the hard scales of the dipoles after subleading emissions.//debug
   */
  void hardScalesSubleading(list<DipoleSplittingInfo> candidates,Energy hardPt);

  /**
   * Return the evolution ordering
   */
  Ptr<DipoleEvolutionOrdering>::tptr evolutionOrdering() const { return theEvolutionOrdering; }

  /**
   * Reshuffle to constituent mass shells
   */
  void constituentReshuffle();

  /**
   * Reshuffle to constituent mass shells
   */
  void decayConstituentReshuffle( PerturbativeProcessPtr decayProc);

  /**
   * Access the generator map
   */
  GeneratorMap& generators() { return theGenerators; }

  /**
   * Access the event record
   */
  DipoleEventRecord& eventRecord() { return theEventRecord; }

  /**
   * Return the event record
   */
  const DipoleEventRecord& eventRecord() const { return theEventRecord; }

  /**
   * Return the splitting kernels.
   */
  const vector<Ptr<DipoleSplittingKernel>::ptr>& splittingKernels() const {
    return kernels;
  }
  
  /**
   * Return the set of offshell parton ids.
   **/
  const set<long>& offShellPartons() { return theColouredOffShellInShower; }
  
  /**
   * Realign the event such as to have the incoming partons along thre
   * beam axes.
   */
  bool realign();

  /**
   * The choice of z boundaries; 0 = restricted, 1 = open, 2 = mixed/other
   */
  virtual int showerPhaseSpaceOption() const {
    return theZBoundaries;
  }

protected:

  /**
   * Perform the cascade.
   */
  void doCascade(unsigned int& emDone,
		 Energy optHardPt = ZERO,
		 Energy optCutoff = ZERO,
		 const bool decay = false);
  /**
   * Set the number of emissions 
   **/
  void setNEmissions(unsigned int n){nEmissions=n;}
  

  /**
   * Get the winning splitting for the
   * given dipole and configuration.
   */
  Energy getWinner(DipoleSplittingInfo& winner,
		   const Dipole& dip,
		   pair<bool,bool> conf,
		   Energy optHardPt = ZERO,
		   Energy optCutoff = ZERO);

  /**
   * Get the winning splitting for the
   * given dipole and configuration.
   */
  Energy getWinner(DipoleSplittingInfo& winner,
		   Energy optHardPt = ZERO,
		   Energy optCutoff = ZERO);

  /**
   * Get the winning splitting for the
   * given dipole and configuration.
   */
  Energy getWinner(SubleadingSplittingInfo& winner,
		   Energy optHardPt = ZERO,
		   Energy optCutoff = ZERO);

  /**
   * Get the winning splitting for the
   * given dipole and configuration.
   */
  Energy getWinner(DipoleSplittingInfo& winner,
		   const DipoleIndex& index,
		   double emitterX, double spectatorX,
		   pair<bool,bool> conf,
		   tPPtr emitter, tPPtr spectator,
		   Energy startScale,
		   Energy optHardPt = ZERO,
		   Energy optCutoff = ZERO);
  
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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


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
   * The splitting kernels to be used.
   */
  vector<Ptr<DipoleSplittingKernel>::ptr> kernels;

  /**
   * The evolution ordering considered
   */
  Ptr<DipoleEvolutionOrdering>::ptr theEvolutionOrdering;

  /**
   * The ConstituentReshuffler to be used
   */
  Ptr<ConstituentReshuffler>::ptr constituentReshuffler;

  /**
   * The intrinsic pt generator to be used.
   */
  Ptr<IntrinsicPtGenerator>::ptr intrinsicPtGenerator;

  /**
   * A global alpha_s to be used for all splitting kernels.
   */
  Ptr<AlphaSBase>::ptr theGlobalAlphaS;

  /**
   * Apply chain ordering to events from matrix
   * element corrections.
   */
  bool chainOrderVetoScales;

  /**
   * Limit the number of emissions.
   * Limit applied if > 0.
   */
  unsigned int nEmissions;

  /**
   * Discard events which did not radiate.
   */
  bool discardNoEmissions;

  /**
   * Perform the first MC@NLO emission only.
   */
  bool firstMCatNLOEmission;

  /**
   * True if powheg style emissions are to be used in the decays
   */
  bool thePowhegDecayEmission;

 /**
   * Switch to record information required for the 
   * nearest neighbour analysis.
   */
  //bool theAnalyseSpinCorrelations;

  /**
   * The realignment scheme
   */
  int realignmentScheme;

  /**
   * Switch on or off subleading Nc corrections
   */
  bool doSubleadingNc;
  
  /**
   * Number of emissions to do subleading Nc corrections to.
   */
  size_t subleadingNcEmissionsLimit;

  /**
   * Current reference weight used for partial unweighting of the 
   * subleading colour shower (updated each emission).
   */
  int currentReferenceWeight;

  /**
   * Integer used to set which method of evolving the density operator
   * to use:
   * 0 - Vijk is Eikonal but there is a cutoff.
   * 1 - Vijk is Eikonal.
   * 2 - Vijk=1 for all i,j,k.
   * 3 - Semi-leading Nc, Vijk=0 for all 
   */
  int densityOperatorEvolution;

  /**
   * Cutoff scale for the invariants (e.g. pEmitter*pEmission) in the 
   * Eikonal dipole kernel, Vijk.
   */
  Energy2 densityOperatorCutoff;

  /**
   * Switch on or off partial unweighting in after each subleading emission.
   */
  bool doPartialUnweightingAtEmission;

  /**
   * Switch on or off partial unweighting in the splitting generator.
   */
  bool doPartialUnweighting;

  /**
   * Reference weight for the partial unweighting.
   */
  double referenceWeight;

  /**
   * Factor changing the acceptance probability for the veto algorithm.
   */
  double cmecReweightFactor;

  /**
   * Scaling factor for the negative colour matrix element corrections.
   */
  double negCMECScaling;

private:

  /**
   * The verbosity level.
   * 0 - print no info
   * 1 - print diagnostic information on setting up
   *     splitting generators etc.
   * 2 - print detailed event information for up to
   *     printEvent events.
   * 3 - print dipole chains after each splitting.
   */
  int verbosity;

  /**
   * See verbosity.
   */
  int printEvent;

private:

  /**
   * The splitting generators indexed by the dipole
   * indices they can work on.
   */
  GeneratorMap theGenerators;

  /**
   * The evnt record used.
   */
  DipoleEventRecord theEventRecord;

  /**
   * The vertex record.
   **/
  DipoleVertexRecord theVertexRecord;

  /**
   * The number of shoer tries so far.
   */
  unsigned int nTries;

  /**
   * Whether or not we did radiate anything
   */
  bool didRadiate;

  /**
   * Whether or not we did realign the event
   */
  bool didRealign;

  /**
   * Vector of candidate splittings containing a vector of the
   * weights, scale and a bool for every step in the reweighted 
   * veto algorithm. The bool is true for an accept step.
   */
  vector<vector<std::tuple<Energy,double,bool> > > theWeightsVector;

  /**
   * Winning candidate index.
   */
  size_t winnerIndex;//debug
  size_t kernelIndex;//debug
  size_t winningKernelIndex;//debug
  vector<Energy> scales;//debug

private:

  /**
   * A freezing value for the renormalization scale
   */
  Energy theRenormalizationScaleFreeze;

  /**
   * A freezing value for the factorization scale
   */
  Energy theFactorizationScaleFreeze;

  /**
   * The matching subtraction, if appropriate
   */
  Ptr<ShowerApproximation>::tptr theShowerApproximation;

  /**
   * True, if sampler should apply compensation
   */
  bool theDoCompensate;

  /**
   * Return the number of accepted points after which the grid should
   * be frozen
   */
  unsigned long theFreezeGrid;

  /**
   * The detuning factor applied to the sampling overestimate kernel
   */
  double theDetuning;

  /**
   * A pointer to the dipole event reweight object
   */
  Ptr<DipoleEventReweight>::ptr theEventReweight;

  /**
   * A pointer to a global dipole splitting reweight
   */
  Ptr<DipoleSplittingReweight>::ptr theSplittingReweight;

  /**
   * True if no warnings have been issued yet
   */
  static bool firstWarn;

  /**
   * The shower starting scale for the last event encountered
   */
  Energy maxPt;

  /**
   * The shower hard scale for the last event encountered
   */
  Energy muPt;
  
  
  /**
   * The merging helper takes care of merging multiple LO and NLO
   * cross sections. Here we need to check if an emission would 
   * radiate in the matrix element region of an other multipicity.
   * If so, the emission is vetoed.
   */
  Ptr<MergerBase>::ptr theMergingHelper;


  /**
   *  PDG codes of the partons which can have an off-shell mass,
   *  this is fast storage for use during running
   */
  set<long> theColouredOffShellInShower;

  /**
   *  PDG codes of the partons which can have an off-shell mass,
   *  this is a vector that is interfaced so they can be changed
   */
  vector<long> theInputColouredOffShellInShower;
  
  /**
   * The choice of z boundaries; 0 = restricted, 1 = open, 2 = mixed/other
   */
  int theZBoundaries;
  
private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DipoleShowerHandler> initDipoleShowerHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleShowerHandler & operator=(const DipoleShowerHandler &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DipoleShowerHandler. */
template <>
struct BaseClassTrait<Herwig::DipoleShowerHandler,1> {
  /** Typedef of the first base class of DipoleShowerHandler. */
  typedef Herwig::ShowerHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DipoleShowerHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DipoleShowerHandler>
  : public ClassTraitsBase<Herwig::DipoleShowerHandler> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DipoleShowerHandler"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DipoleShowerHandler is implemented. It may also include several, space-separated,
   * libraries if the class DipoleShowerHandler depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DipoleShowerHandler_H */
